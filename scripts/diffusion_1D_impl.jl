using Plots, Printf, LinearAlgebra

@views function diffusion_1D()
    # Physics
    lx    = 10.0       # domain size
    D     = 1.0        # diffusion coefficient
    ttot  = 0.6        # total simulation time
    dt    = 0.002      # physical time step
    # Numerics
    nx    = 128        # numerical grid resolution
    tol   = 1e-4       # tolerance
    itMax = 1e5        # max number of iterations
    # Derived numerics
    dx    = lx/nx      # grid size
    dtau  = (1.0/(dx^2/D/2.1) + 1.0/dt)^-1 # iterative "timestep"
    xc    = LinRange(dx/2, lx-dx/2, nx)
    # Array allocation
    qH    = zeros(nx-1)
    dHdtau  = zeros(nx-2)
    # Initial condition
    H0    = exp.(-(xc.-lx/2).^2)
    Hold  = copy(H0)
    H     = copy(H0)
    t = 0.0; it = 1; err = 2*tol
    # Physical time loop
    while t<ttot
        # Picard-type iteration
        while err>tol && it<itMax
            qH         .= -D*diff(H)/dx            # flux
            dHdtau     .= -(H[2:end-1] - Hold[2:end-1])/dt - diff(qH)/dx  # rate of change
            H[2:end-1] .= H[2:end-1] + dtau*dHdtau  # update rule
            it += 1; err = norm(dHdtau)/length(dHdtau)
        end
        err = 2*tol; t += dt
        Hold .= H
    end
    # Analytic solution
    Hana     = 1/sqrt(4*(ttot+1/4)) * exp.(-(xc.-lx/2).^2 /(4*(ttot+1/4)))
    @show norm(H.-Hana, Inf), norm(H.-Hana, Inf)<0.00084

    @printf("Total time = %1.2f, it tot = %d \n", round(dt, sigdigits=2), it)
    # Visualize
    plot(xc,0* H0, linewidth=3); display(plot!(xc, H.-Hana, legend=false, framestyle=:box, linewidth=3, xlabel="lx", ylabel="H", title="implicit diffusion (niter=$it)"))

    return
end

diffusion_1D()
