using Plots, Printf

@views function iceflow()
    # physics
    s2y    = 3600*24*365.25  # seconds to years
    lx, ly = 30e3, 30e3
    rho_i  = 910.0           # ice density
    g      = 9.81            # gravity acceleration
    npow   = 3.0             # Glen's power law exponent 
    a0     = 1.5e-24         # Glen's law enhancement term
    # numerics
    nx, ny = 100, 100
    nt     = 1e5
    nout   = 200
    tolnl  = 1e-8
    epsi   = 1e-2
    damp   = 0.92
    # derived physics
    a      = 2.0*a0/(npow+2)*(rho_i*g)^npow*s2y
    # derived numerics
    dx, dy = lx/nx, ly/ny
    xc, yc = LinRange(dx/2, lx-dx/2, nx), LinRange(dy/2, ly-dy/2, ny)
    (Xc,Yc)= ([x for x=xc,y=yc], [y for x=xc,y=yc])
    cfl    = max(dx^2,dy^2)/4.1/2.0
    # array initialisation
    S      = zeros(nx  , ny  )
    B      = zeros(nx  , ny  )
    H      =  ones(nx  , ny  )
    M      = zeros(nx  , ny  )
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
    # initial condition (Jarosch's 2013 benchmark)
    xm     = 20e3
    xmB    = 7e3
    M     .= (((npow.*2.0./xm.^(2*npow-1)).*Xc.^(npow-1)).*abs.(xm.-Xc).^(npow-1)).*(xm.-2.0.*Xc)
    M[Xc.>xm]  .= 0.0
    B[Xc.<xmB] .= 500
    # smoothing (Mahaffy, 1976)
    B[2:end-1,2:end-1] .= B[2:end-1,2:end-1] .+ 1.0./4.1.*(diff(diff(B[:,2:end-1], dims=1), dims=1) .+ diff(diff(B[2:end-1,:], dims=2), dims=2))
    S     .= B .+ H
    # time loop
    for it = 1:nt
        Err   .= H
        # compute diffusivity
        dSdx  .= diff(S, dims=1)/dx
        dSdy  .= diff(S, dims=2)/dy
        gradS .= sqrt.(av_ya(dSdx).^2 .+ av_xa(dSdy).^2)
        D     .= a*av(H).^(npow+2) .* gradS.^(npow-1)
        # compute flux
        qHx   .= .-av_ya(D).*diff(S[:,2:end-1], dims=1)/dx
        qHy   .= .-av_xa(D).*diff(S[2:end-1,:], dims=2)/dy
        # update ice thickness
        dt    .= 0.5*min.(1.0, cfl./(epsi .+ av(D)))
        ResH  .= .-(diff(qHx, dims=1)/dx .+ diff(qHy, dims=2)/dy) .+ inn(M)
        dHdt  .= dHdt.*damp .+ ResH
        H[2:end-1,2:end-1] .= max.(0.0, inn(H) .+ dt.*dHdt)
        # update surface
        S     .= B .+ H
        if mod(it, nout)==0
            # error check
            Err  .= Err .- H
            err   = (sum(abs.(Err))./nx./ny)
            @printf("it = %d, error = %1.2e \n", it, err)
            # stop criterion
            if (err<tolnl) break end
        end
    end
    # compute velocities
    Vx .= -D./(av(H) .+ epsi).*av_ya(dSdx)
    Vy .= -D./(av(H) .+ epsi).*av_xa(dSdy)
    # visualisation
    xc, yc = xc./1e3, yc./1e3
    xv, yv = 0.5*(xc[1:end-1].+xc[2:end]), 0.5*(yc[1:end-1].+yc[2:end])
    p1 = heatmap(xc, yc, S' , aspect_ratio=1, xlims=(xc[1], xc[end]), ylims=(yc[1], yc[end]), c=:davos, title="S")
    p2 = heatmap(xc, yc, H' , aspect_ratio=1, xlims=(xc[1], xc[end]), ylims=(yc[1], yc[end]), c=:davos, title="H")
    p3 = heatmap(xv, yv, Vx', aspect_ratio=1, xlims=(xv[1], xv[end]), ylims=(yv[1], yv[end]), c=:davos, title="Vx")
    p4 = heatmap(xv, yv, Vy', aspect_ratio=1, xlims=(xv[1], xv[end]), ylims=(yv[1], yv[end]), c=:davos, title="Vy")
    display(plot( p1, p2, p3, p4))
    return
end

iceflow()
