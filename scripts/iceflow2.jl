using NCDatasets, Plots, Printf

@views av(A)    = 0.25*(A[1:end-1,1:end-1].+A[2:end,1:end-1].+A[1:end-1,2:end].+A[2:end,2:end])
@views av_xa(A) = 0.5.*(A[1:end-1,:].+A[2:end,:])
@views av_ya(A) = 0.5.*(A[:,1:end-1].+A[:,2:end])
@views inn(A)   = A[2:end-1,2:end-1]

@views function get_data(url, downscale=30)
    # Download BedMachine v3
    !ispath("../data") && mkdir("../data")
    !ispath("../output") && mkdir("../output")
    filename = joinpath("../data", splitdir(url)[2])
    print("Loading NC file ... ")
    ds = NCDataset(filename)
    # mask: 0=ocean, 1=ice-free land, 2=grounded-ice, 3=floating-ice, 4=non-Greenland land
    x, y   = ds["x"][1:downscale:end], ds["y"][1:downscale:end]
    x, y   = x[1]:x[2]-x[1]:x[end], y[1]:y[2]-y[1]:y[end]
    dx, dy = x[2]-x[1], y[2]-y[1]
    mask_     = ds["mask"][1:downscale:end, 1:downscale:end]
    # mask      = (mask_.==0) .| (mask_.==1) .| (mask_.==2) .| (mask_.==3)
    mask      = (mask_.==1) .| (mask_.==2)
    surface   = convert(Matrix{Float64}, ds["surface"][1:downscale:end, 1:downscale:end])
    thickness = convert(Matrix{Float64}, ds["thickness"][1:downscale:end, 1:downscale:end])
    bed       = convert(Matrix{Float64}, ds["bed"][1:downscale:end, 1:downscale:end])
    println("done")
    # Pad the data to make it GPU compatible with good performance # DEBUG: do this until ParallelStencil is fixed: see https://github.com/samo-lin/ParallelStencil.jl/issues/42
    print("Padding the data for GPUs optimal perf ... ")
    olx, oly = rem(length(x), 16), rem(length(y), 16)
    x1=1; xE=0; if (olx!=0) x1 = Int(ceil(olx/2)); xE = olx-x1+1; end
    y1=1; yE=0; if (oly!=0) y1 = Int(ceil(oly/2)); yE = oly-y1+1; end
    x, y   = x[x1:end-xE], y[y1:end-yE]
    dx, dy = x[2]-x[1], y[2]-y[1]
    mask      =      mask[x1:end-xE, y1:end-yE]
    surface   =   surface[x1:end-xE, y1:end-yE]
    thickness = thickness[x1:end-xE, y1:end-yE]
    bed       =       bed[x1:end-xE, y1:end-yE]
    println("done")
    println("Grid size: nx=$(length(x)), ny=$(length(y))")
    return mask, surface, thickness, bed, x, y, dx, dy
end

@views function iceflow(dx, dy, Zbed, Hice, Mask=zero(Zbed); do_visu=false)
    print("Starting ice flow model ... ")
    # physics
    s2y    = 3600*24*365.25  # seconds to years
    rho_i  = 910.0           # ice density
    g      = 9.81            # gravity acceleration
    npow   = 3.0             # Glen's power law exponent 
    a0     = 1.5e-24         # Glen's law enhancement term
    b_max  = 0.13
    # b = min.(((1.3517 - 0.014158*82)/100*0.91).*(Z.-800),b_max); # doi: 10.1017/jog.2016.75
    # numerics
    nx, ny = size(Zbed,1), size(Zbed,2)
    @assert (nx, ny) == size(Zbed) == size(Hice) == size(Mask) "Sizes don't match"
    nt     = 1e5
    nout   = 200
    tolnl  = 1e-6
    epsi   = 1e-4
    damp   = 0.85
    ns     = 2
    # derived physics
    a      = 2.0*a0/(npow+2)*(rho_i*g)^npow*s2y
    lx, ly = nx*dx, ny*dy
    # derived numerics
    xc, yc = LinRange(dx/2, lx-dx/2, nx), LinRange(dy/2, ly-dy/2, ny)
    xv, yv = 0.5*(xc[1:end-1].+xc[2:end]), 0.5*(yc[1:end-1].+yc[2:end])
    (Xc,Yc)= ([x for x=xc,y=yc], [y for x=xc,y=yc])
    cfl    = max(dx^2,dy^2)/4.1
    dtsc   = 1.0/3.0
    # array initialisation
    Err    = zeros(nx  , ny  )
    dSdx   = zeros(nx-1, ny  )
    dSdy   = zeros(nx  , ny-1)
    gradS  = zeros(nx-1, ny-1)
    D      = zeros(nx-1, ny-1)
    qHx    = zeros(nx-1, ny-2)
    qHy    = zeros(nx-2, ny-1)
    dt     = zeros(nx-2, ny-2)
    ResH   = zeros(nx-2, ny-2)
    dHdt   = zeros(nx-2, ny-2)
    Vx     = zeros(nx-1, ny-1)
    Vy     = zeros(nx-1, ny-1)
    M      = zeros(nx  , ny  )
    # initial condition
    S      = zeros(nx  , ny  )
    B      = Zbed #zeros(nx  , ny  )
    H      = Hice # ones(nx  , ny  )
    Yc2    = Yc .- minimum(Yc); Yc2 .= Yc2./maximum(Yc2)
    grad_b = (1.3517 .- 0.014158.*(60.0.+Yc2*20.0))./100.0.*0.91
    z_ELA  = 1500.0 .- Yc2*400.0
    # display(heatmap(xc./1e3, reverse(yc)./1e3, reverse(grad_b, dims=2)', c=:davos, aspect_ratio=1, xlims=(xc[1], xc[end])./1e3, ylims=(yc[end], yc[1])./1e3, framestyle=:box, title="Surface"))
    # error("stop")
    if do_visu
        H_v = fill(NaN, nx, ny)
        S_v = fill(NaN, nx, ny)
        M_v = fill(NaN, nx, ny)
        V_v = fill(NaN, nx-2, ny-2)
    end
    # smoothing (Mahaffy, 1976)
    for is=1:ns
        B[2:end-1,2:end-1] .= B[2:end-1,2:end-1] .+ 1.0./4.1.*(diff(diff(B[:,2:end-1], dims=1), dims=1) .+ diff(diff(B[2:end-1,:], dims=2), dims=2))
        H[2:end-1,2:end-1] .= H[2:end-1,2:end-1] .+ 1.0./4.1.*(diff(diff(H[:,2:end-1], dims=1), dims=1) .+ diff(diff(H[2:end-1,:], dims=2), dims=2))
    end
    S     .= B .+ H
    println("time loop:")
    # time loop
    for it = 1:nt
        Err   .= H
        # mass balance
        M     .= min.(grad_b.*(S.-z_ELA), b_max)
        # compute diffusivity
        dSdx  .= diff(S, dims=1)/dx
        dSdy  .= diff(S, dims=2)/dy
        gradS .= sqrt.(av_ya(dSdx).^2 .+ av_xa(dSdy).^2)
        D     .= a*av(H).^(npow+2) .* gradS.^(npow-1)
        # compute flux
        qHx   .= .-av_ya(D).*diff(S[:,2:end-1], dims=1)/dx
        qHy   .= .-av_xa(D).*diff(S[2:end-1,:], dims=2)/dy
        # update ice thickness
        # dt    .= dtsc*min.(1.0, cfl./(epsi .+ av(D)))
        dt    .= dtsc*min.(10.0, cfl./(epsi .+ av(D)))
        ResH  .= .-(diff(qHx, dims=1)/dx .+ diff(qHy, dims=2)/dy) .+ inn(M)
        dHdt  .= dHdt.*damp .+ ResH
        H[2:end-1,2:end-1] .= max.(0.0, inn(H) .+ dt.*dHdt)
        # apply mask
        H[Mask.==0] .= 0.0
        # update surface
        S     .= B .+ H
        if mod(it, nout)==0
            # error check
            Err .= Err .- H
            err  = (sum(abs.(Err))./nx./ny)
            @printf("it = %d, error = %1.2e \n", it, err)
            # stop criterion
            if (err<tolnl) break end
            p1 = heatmap(xc./1e3, reverse(yc)./1e3, reverse(H, dims=2)', c=:davos, aspect_ratio=1, xlims=(xc[1], xc[end])./1e3, ylims=(yc[end], yc[1])./1e3, framestyle=:box, title="Ice thickness")
            p2 = heatmap(xc[2:end-1]./1e3, reverse(yc[2:end-1])./1e3, reverse(dt, dims=2)', c=:davos, aspect_ratio=1, xlims=(xc[2], xc[end-1])./1e3, ylims=(yc[end-1], yc[2])./1e3, framestyle=:box, title="Surface")
            display(plot(p1, p2))
        end
    end
    # compute velocities
    Vx .= -D./(av(H) .+ epsi).*av_ya(dSdx)
    Vy .= -D./(av(H) .+ epsi).*av_xa(dSdy)
    # visualisation
    if do_visu
        H_v.=H; H_v[H.==0].=NaN
        S_v.=S; S_v[Mask.==0].=NaN
        M_v.=M; M_v[Mask.==0].=NaN
        V_v.=sqrt.(av(Vx).^2 .+ av(Vy).^2); V_v[inn(H).==0].=NaN
        p1 = heatmap(xc./1e3, reverse(yc)./1e3, reverse(S_v, dims=2)', c=:davos, aspect_ratio=1, xlims=(xc[1], xc[end])./1e3, ylims=(yc[end], yc[1])./1e3, framestyle=:box, title="Surface")
        p2 = heatmap(xc./1e3, reverse(yc)./1e3, reverse(H_v, dims=2)', c=:davos, aspect_ratio=1, xlims=(xc[1], xc[end])./1e3, ylims=(yc[end], yc[1])./1e3, framestyle=:box, title="Ice thickness")
        p3 = heatmap(xc[2:end-1]./1e3, reverse(yc[2:end-1])./1e3, reverse(V_v, dims=2)', c=:hot, aspect_ratio=1, xlims=(xc[2], xc[end-1])./1e3, ylims=(yc[end-1], yc[2])./1e3, clims=(0.0, 1000.0), framestyle=:box, title="||vel||")
        p4 = heatmap(xc./1e3, reverse(yc)./1e3, reverse(M_v, dims=2)', c=:davos, aspect_ratio=1, xlims=(xc[1], xc[end])./1e3, ylims=(yc[end], yc[1])./1e3, framestyle=:box, title="Mass Bal.")
        display(plot(p1, p2, p3, p4))
    end
    return H, S, M, Vx, Vy
end

# ------------------------------------------------------------------------------

Mask, Surf, Hice, Zbed, xc, yc, dx, dy = get_data("BedMachineGreenland-2017-09-20.nc", 60)
# xc, yc = xc./1e3, yc./1e3
# xv, yv = 0.5*(xc[1:end-1].+xc[2:end]), 0.5*(yc[1:end-1].+yc[2:end])
# p1 = heatmap(xc,reverse(yc),reverse(Surf, dims=2)', c=:davos, title="Surf")
# p2 = heatmap(xc,reverse(yc),reverse(Hice, dims=2)', c=:davos, title="Hice")
# p3 = heatmap(xc,reverse(yc),reverse(Zbed, dims=2)', c=:davos, title="Zbed")
# p4 = heatmap(xc,reverse(yc),reverse(Mask, dims=2)', c=:davos, title="Mask")
# display(plot(p1, p2, p3, p4))

H, S, M, Vx, Vy = iceflow(dx, dy, Zbed, Surf.-Zbed, Mask; do_visu=true)

H_diff = Hice-H; H_diff[Mask.==0] .= NaN
Hice[Mask.==0] .= NaN
H[Mask.==0] .= NaN

p1 = heatmap(xc./1e3, reverse(yc)./1e3, reverse(Hice, dims=2)', c=:davos, aspect_ratio=1, xlims=(xc[1], xc[end])./1e3, ylims=(yc[end], yc[1])./1e3, framestyle=:box, title="Hdata")
p2 = heatmap(xc./1e3, reverse(yc)./1e3, reverse(H, dims=2)', c=:davos, aspect_ratio=1, xlims=(xc[1], xc[end])./1e3, ylims=(yc[end], yc[1])./1e3, framestyle=:box, title="Hmodel")
p3 = heatmap(xc./1e3, reverse(yc)./1e3, reverse(H_diff, dims=2)', c=:davos, aspect_ratio=1, xlims=(xc[1], xc[end])./1e3, ylims=(yc[end], yc[1])./1e3, framestyle=:box, title="Hdata-Hmodel")
display(plot(p1, p2, p3, layout=(1, 3)))

println("... done.")
