const USE_GPU = false
using ParallelStencil
using ParallelStencil.FiniteDifferences1D
@static if USE_GPU
    @init_parallel_stencil(CUDA, Float64, 1)
else
    @init_parallel_stencil(Threads, Float64, 1)
end
using Plots, Printf

@parallel function compute_flux!(qH, H, D, dx)
    @all(qH) = -D*@d(H)/dx
    return
end

@parallel function compute_rate!(dHdt, H, Hold, qH, dt, damp, dx)
    @all(dHdt) = -(@all(H) - @all(Hold))/dt -@d(qH)/dx + damp*@all(dHdt)
    return
end

@parallel function compute_update!(H, dHdt, dtau)
    @inn(H) = @inn(H) + dtau*@all(dHdt)
    return
end

@views function diffusion_1D()
    # Physics
    lx   = 10.0       # domain size
    D    = 1.0        # diffusion coefficient
    dt   = 0.6        # physical time step (ttot)
    # Numerics
    nx   = 128        # numerical grid resolution
    epsi = 1e-6       # tolerance
    damp = 0.86       # damping
    # Derived numerics
    dx   = lx/nx      # grid size
    dtau = (1.0/(dx^2/D/2.1) + 1.0/dt)^-1 # iterative timestep
    xc   = LinRange(dx/2, lx-dx/2, nx)
    # Array allocation
    qH   = @zeros(nx-1)
    dHdt = @zeros(nx-2)
    # Initial condition
    H    = Data.Array( exp.(.-(xc.-lx./2.0).^2) )
    Hold = @ones(nx).*H
    H0   = @ones(nx).*H
    it = 1; err = 2*epsi
    # Time loop
    while err>epsi
        @parallel compute_flux!(qH, H, D, dx)
        @parallel compute_rate!(dHdt, H, Hold, qH, dt, damp, dx)
        @parallel compute_update!(H, dHdt, dtau)
        it += 1; err = maximum(dHdt)
    end
    @printf("Total time = %1.2f, it tot = %d \n", round(dt, sigdigits=2), it)
    # Visualise
    plot(xc, H0, linewidth=3); display(plot!(xc, H, legend=false, framestyle=:box, linewidth=3, xlabel="lx", ylabel="H", title="damped diffusion (niter=$it)"))
    return
end

diffusion_1D()
