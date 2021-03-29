using Plots, Printf

@views function diffusion_1D()
    # Physics
    lx    = 10.0        # domain size
    D     = 1.0         # diffusion coefficient
    ttot  = 0.6         # total simulation time
    # Numerics
    nx    = 128         # numerical grid resolution
    # Derived numerics
    dx    = lx/nx       # gird size
    dt    = dx^2/D/2.1  # time step
    xc    = LinRange(dx/2, lx-dx/2, nx)
    # Array allocation
    qH    = zeros(nx-1)
    dHdt  = zeros(nx-2)
    # Initial condition
    H     = exp.(.-(xc.-lx./2.0).^2)
    H0    = copy(H)
    t = 0.0; it = 1
    # Time loop
    while t<ttot
        qH         .= -D*diff(H)/dx         # flux
        dHdt       .=  -diff(qH)/dx         # rate of change
        H[2:end-1] .= H[2:end-1] .+ dt*dHdt # update rule
        t += dt; it += 1
    end
    @printf("Total time = %1.2f, it tot = %d \n", round(ttot, sigdigits=2), it)
    # Visualise
    plot(xc, H0, linewidth=3); display(plot!(xc, H, legend=false, framestyle=:box, linewidth=3, xlabel="lx", ylabel="H", title="explicit diffusion (nt=$it)"))
    return
end

diffusion_1D()
