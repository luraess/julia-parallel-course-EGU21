using Plots, Printf

@views function diffusion_1D()
    # Physics
    lx   = 10.0       # domain size
    D    = 1.0        # diffusion coefficient
    dt   = 0.6        # physical time step (ttot)
    # Numerics
    nx   = 128        # numerical grid resolution
    epsi = 1e-6       # tolerance
    # Derived numerics
    dx   = lx/nx      # grid size
    dtau = (1.0/(dx^2/D/2.1) + 1.0/dt)^-1 # iterative timestep
    xc   = LinRange(dx/2, lx-dx/2, nx)
    # Array allocation
    qH   = zeros(nx-1)
    dHdt = zeros(nx-2)
    # Initial condition
    H    = exp.(.-(xc.-lx./2.0).^2)
    Hold = copy(H)
    H0   = copy(H)
    it = 1; err = 2*epsi
    # Time loop
    while err>epsi
        qH         .= -D*diff(H)/dx            # flux
        dHdt       .= -(H[2:end-1].-Hold[2:end-1])/dt .-diff(qH)/dx  # rate of change
        H[2:end-1] .= H[2:end-1] .+ dtau*dHdt  # update rule
        it += 1; err = maximum(dHdt)
    end
    @printf("Total time = %1.2f, it tot = %d \n", round(dt, sigdigits=2), it)
    # Visualise
    plot(xc, H0, linewidth=3); display(plot!(xc, H, legend=false, framestyle=:box, linewidth=3, xlabel="lx", ylabel="H", title="implicit diffusion (niter=$it)"))
    return
end

diffusion_1D()
