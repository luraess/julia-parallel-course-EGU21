const USE_GPU = false
using ParallelStencil
using ParallelStencil.FiniteDifferences2D
@static if USE_GPU
    @init_parallel_stencil(CUDA, Float64, 2)
    macro pow(args...)  esc(:(CUDA.pow($(args...)))) end
else
    @init_parallel_stencil(Threads, Float64, 2)
    pow(x,y) = x^y
    macro pow(args...)  esc(:(pow($(args...)))) end
end
using JLD, Plots, Printf, LinearAlgebra

# CPU functions
@views av(A)  = 0.25*(A[1:end-1,1:end-1].+A[2:end,1:end-1].+A[1:end-1,2:end].+A[2:end,2:end])
@views inn(A) = A[2:end-1,2:end-1]

@views function apply_smoothing!(A, B, ns)
    print("Applying some smoothing ... ")
    A2 = @zeros(size(A)); A2 .= A
    B2 = @zeros(size(B)); B2 .= B
    for is=1:ns
        @parallel smooth!(A2, A)
        @parallel smooth!(B2, B)
        @parallel (1:size(A2,2)) bc_x!(A2)
        @parallel (1:size(A2,1)) bc_y!(A2)
        @parallel (1:size(B2,2)) bc_x!(B2)
        @parallel (1:size(B2,1)) bc_y!(B2)
        A, A2 = A2, A
        B, B2 = B2, B
    end
    println("done.")
    return
end

# GPU function
@parallel function smooth!(A2, A)
    @inn(A2) = @inn(A) + 1.0./4.1.*(@d2_xi(A) + @d2_yi(A))
    return
end

@parallel_indices (iy) function bc_x!(A::Data.Array)
    A[1  , iy] = A[2    , iy]
    A[end, iy] = A[end-1, iy]
    return
end

@parallel_indices (ix) function bc_y!(A::Data.Array)
    A[ix, 1  ] = A[ix, 2    ]
    A[ix, end] = A[ix, end-1]
    return
end

@parallel function compute_Err1!(Err, H)
    @all(Err) = @all(H)
    return
end

@parallel function compute_Err2!(Err, H)
    @all(Err) = @all(Err) - @all(H)
    return
end

@parallel function assign!(Aold, A)
    @all(Aold) = @all(A)
    return
end

@parallel function compute_M_dS!(M, dSdx, dSdy, S, z_ELA, grad_b, b_max, dELA, dx, dy)
    @all(M)    = min(@all(grad_b)*(@all(S) - (@all(z_ELA) + dELA)), b_max)
    @all(dSdx) = @d_xa(S)/dx
    @all(dSdy) = @d_ya(S)/dy
    return
end

@parallel function compute_D!(D, H, dSdx, dSdy, a, npow)
    @all(D) = a*@pow(@av(H), (npow+2)) * (@av_ya(dSdx)*@av_ya(dSdx) + @av_xa(dSdy)*@av_xa(dSdy))
    return
end

@parallel function compute_flux!(qHx, qHy, D, S, dt, dx, dy)
    @all(qHx)  = -@av_ya(D)*@d_xi(S)/dx
    @all(qHy)  = -@av_xa(D)*@d_yi(S)/dy
    
    return
end

@parallel function compute_dHdt!(dtau, ResH, dHdt, D, H, Hold, qHx, qHy, M, dtausc, cfl, epsi, dt, damp, dx, dy)
    @all(dtau) = dtausc*min(10.0, 1.0/(1.0/dt + 1.0/(cfl/(epsi + @av(D)))))
    @all(ResH) = -(@inn(H) - @inn(Hold))/dt -(@d_xa(qHx)/dx + @d_ya(qHy)/dy) + @inn(M)
    @all(dHdt) = @all(dHdt)*damp + @all(ResH)
    return
end

@parallel function compute_H!(H, dHdt, dtau)
    @inn(H) = max(0.0, @inn(H) + @all(dtau)*@all(dHdt) )
    return
end

@parallel_indices (ix,iy) function compute_Mask_S!(H, S, B, Mask)
    if (ix<=size(H,1) && iy<=size(H,2)) if (Mask[ix,iy]==0) H[ix,iy] = 0.0 end end
    if (ix<=size(S,1) && iy<=size(S,2)) S[ix,iy] = B[ix,iy] + H[ix,iy] end    
    return
end

@parallel function compute_Vel!(Vx, Vy, D, H, dSdx, dSdy, epsi)
    @all(Vx) = -@all(D)/(@av(H) + epsi)*@av_ya(dSdx)
    @all(Vy) = -@all(D)/(@av(H) + epsi)*@av_xa(dSdy)
    return
end

@views function iceflow(dx, dy, Zbed, Hice, Mask; do_visu=false)
    println("Initialising ice flow model ... ")
    # physics
    s2y      = 3600*24*365.25  # seconds to years
    rho_i    = 910.0           # ice density
    g        = 9.81            # gravity acceleration
    npow     = 3.0             # Glen's power law exponent
    a0       = 1.5e-24         # Glen's law enhancement term
    b_max    = 0.15            # max. Mass balance rate
    dt0      = 50.0            # physical timestep in years
    ttot     = 50*dt0          # total time in years
    dELA     = 1               # ELA shift in [m/yr] (~0.6°C / 100m) -> +1.2°C == 200m 
    # numerics
    nx, ny   = size(Zbed,1), size(Zbed,2) # numerical grid resolution
    @assert (nx, ny) == size(Zbed) == size(Hice) == size(Mask) "Sizes don't match"
    itMax    = 1e5             # number of iteration (max)
    nout     = 200             # error check frequency
    tolnl    = 1e-6            # nonlinear tolerance
    epsi     = 1e-4            # small number
    damp     = 0.85            # convergence accelerator
    dtausc   = 1.0/3.0         # iterative dtau scaling
    # derived physics
    a        = 2.0*a0/(npow+2)*(rho_i*g)^npow*s2y
    lx, ly   = nx*dx, ny*dy
    # derived numerics
    xc, yc   = LinRange(dx/2, lx-dx/2, nx), LinRange(dy/2, ly-dy/2, ny)
    xv, yv   = 0.5*(xc[1:end-1].+xc[2:end]), 0.5*(yc[1:end-1].+yc[2:end])
    (Xc,Yc)  = ([x for x=xc,y=yc], [y for x=xc,y=yc])
    cfl      = max(dx^2,dy^2)/4.1
    # array initialisation
    Hold     = @zeros(nx  , ny  )
    Err      = @zeros(nx  , ny  )
    dSdx     = @zeros(nx-1, ny  )
    dSdy     = @zeros(nx  , ny-1)
    D        = @zeros(nx-1, ny-1)
    qHx      = @zeros(nx-1, ny-2)
    qHy      = @zeros(nx-2, ny-1)
    dtau     = @zeros(nx-2, ny-2)
    ResH     = @zeros(nx-2, ny-2)
    dHdt     = @zeros(nx-2, ny-2)
    Vx       = @zeros(nx-1, ny-1)
    Vy       = @zeros(nx-1, ny-1)
    M        = @zeros(nx  , ny  )
    B        = @zeros(nx  , ny  )
    H        = @zeros(nx  , ny  )
    # initial condition
    S        = @zeros(nx  , ny  )
    B       .= Zbed
    H       .= Hice
    Mask     = Data.Array(Mask)
    Yc2      = Yc .- minimum(Yc); Yc2 .= Yc2./maximum(Yc2)
    grad_b   = Data.Array((1.3517 .- 0.014158.*(60.0.+Yc2*20.0))./100.0.*0.91) # Mass Bal. gradient, from doi: 10.1017/jog.2016.75
    z_ELA    = Data.Array(1300.0 .- Yc2*300.0)                                 # Educated guess for ELA altitude
    S       .= B .+ H
    ELA      = 0.0
    if do_visu
        FS = 7
        H_v = fill(NaN, nx, ny)
        S_v = fill(NaN, nx, ny)
        M_v = fill(NaN, nx, ny)
        V_v = fill(NaN, nx-2, ny-2)
        Ma  = Array(Mask)
    end
    println(" starting time loop:")
    t = 0.0; it = 0
    while t<ttot
        println(" > step $(it)")
        if it==0
            dt = 1e50
        else
            dt  = dt0
            ELA = ELA + dt*dELA
        end
        @parallel assign!(Hold, H)
        # iteration loop
        iter = 1; err = 2*tolnl
        while err>tolnl && iter<itMax
            @parallel compute_Err1!(Err, H) 
            @parallel compute_M_dS!(M, dSdx, dSdy, S, z_ELA, grad_b, b_max, ELA, dx, dy)
            @parallel compute_D!(D, H, dSdx, dSdy, a, npow)
            @parallel compute_flux!(qHx, qHy, D, S, dt, dx, dy)
            @parallel compute_dHdt!(dtau, ResH, dHdt, D, H, Hold, qHx, qHy, M, dtausc, cfl, epsi, dt, damp, dx, dy)
            @parallel compute_H!(H, dHdt, dtau)
            @parallel compute_Mask_S!(H, S, B, Mask)
            # error check
            if mod(iter, nout)==0
                @parallel compute_Err2!(Err, H)
                err = norm(Err)/length(Err)
                @printf(" iter = %d, error = %1.2e \n", iter, err)
                if isnan(err) error("NaNs") end # safeguard
            end
            iter += 1
        end
        @parallel compute_Vel!(Vx, Vy, D, H, dSdx, dSdy, epsi)
        if (it>0)  t += dt  end
        # visualisation
        H_v .= Array(H); H_v[Ma.==0] .= NaN
        S_v .= Array(S); S_v[Ma.==0] .= NaN
        M_v .= Array(M); M_v[Ma.==0] .= NaN
        V_v .= sqrt.(av(Array(Vx)).^2 .+ av(Array(Vy)).^2); V_v[inn(Array(H)).==0] .= NaN
        p1 = heatmap(xc./1e3, reverse(yc)./1e3, reverse(S_v, dims=2)', c=:davos, aspect_ratio=1, xlims=(xc[1], xc[end])./1e3, ylims=(yc[end], yc[1])./1e3, yaxis=font(FS, "Courier"), ticks=nothing, framestyle=:box, title="Surface elev. [m]", titlefontsize=FS, titlefont="Courier")
        p2 = heatmap(xc./1e3, reverse(yc)./1e3, reverse(H_v, dims=2)', c=:davos, aspect_ratio=1, xlims=(xc[1], xc[end])./1e3, ylims=(yc[end], yc[1])./1e3, yaxis=font(FS, "Courier"), ticks=nothing, framestyle=:box, title="$(round(t,sigdigits=2)) yrs, Ice thick. [m]", titlefontsize=FS, titlefont="Courier")
        p3 = heatmap(xc[2:end-1]./1e3, reverse(yc[2:end-1])./1e3, reverse(log10.(V_v), dims=2)', c=:batlow, aspect_ratio=1, xlims=(xc[2], xc[end-1])./1e3, ylims=(yc[end-1], yc[2])./1e3, clims=(0.1, 2.0), yaxis=font(FS, "Courier"), ticks=nothing, framestyle=:box, title="log10(vel) [m/yr]", titlefontsize=FS, titlefont="Courier")
        p4 = heatmap(xc./1e3, reverse(yc)./1e3, reverse(M_v, dims=2)', c=:devon, aspect_ratio=1, xlims=(xc[1], xc[end])./1e3, ylims=(yc[end], yc[1])./1e3, yaxis=font(FS, "Courier"), ticks=nothing, framestyle=:box, title="Mass Bal. rate [m/yr]", titlefontsize=FS, titlefont="Courier")
        # display(plot(p1, p2, p3, p4, size=(400,400)))
        plot(p1, p2, p3, p4, size=(400,400), dpi=200); frame(anim) #background_color=:transparent, foreground_color=:white
        # savefig("../output_evo/iceflow_out1_xpu_evo_$(nx)x$(ny)_$(it).png")
        it += 1
    end
    return Array(H), Array(S), Array(M), Array(Vx), Array(Vy)
end
# ------------------------------------------------------------------------------
# load the data
print("Loading the data ... ")
# data = load("../data/BedMachineGreenland_96_184_ds100.jld")
data = load("../data/BedMachineGreenland_160_304_ds60.jld")
Hice, Mask, Zbed = Data.Array(data["Hice"]), Data.Array(data["Mask"]), Data.Array(data["Zbed"])
xc, yc, dx, dy   = data["xc"], data["yc"], data["dx"], data["dy"]
println("done.")

# apply some smoothing
ns = 2
apply_smoothing!(Hice, Zbed, ns)

# handle output
do_visu = true
do_save = true

nx, ny = size(Hice)
if do_visu
    ENV["GKSwstype"]="nul"
    !ispath("../output_evo") && mkdir("../output_evo")
    !ispath("../output_evo/gif_$(nx)x$(ny)") && mkdir("../output_evo/gif_$(nx)x$(ny)"); loadpath = "../output_evo/gif_$(nx)x$(ny)"; anim = Animation(loadpath,String[])
    println("Animation directory: $(anim.dir)")
end

# run the SIA flow model
H, S, M, Vx, Vy = iceflow(dx, dy, Zbed, Hice, Mask; do_visu)

# output and save
if do_visu gif(anim, "../output_evo/iceflow_evo_$(nx)x$(ny).gif", fps = 5) end

if do_save
    save("../output_evo/iceflow_xpu_evo_$(nx)x$(ny).jld", "Hice", convert(Matrix{Float32}, Hice),
                                                          "Mask", convert(Matrix{Float32}, Mask),
                                                          "H"   , convert(Matrix{Float32}, H),
                                                          "S"   , convert(Matrix{Float32}, S),
                                                          "M"   , convert(Matrix{Float32}, M),
                                                          "Vx"  , convert(Matrix{Float32}, Vx),
                                                          "Vy"  , convert(Matrix{Float32}, Vy),
                                                          "xc", xc, "yc", yc)
end

println("... done.")
