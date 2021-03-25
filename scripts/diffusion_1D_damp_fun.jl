using Plots, Printf

function compute_flux!(qH, H, D, dx, nx)
    Threads.@threads for ix=1:nx
        if (ix<=nx-1)  qH[ix] = -D*(H[ix+1]-H[ix])/dx  end
    end
    return
end

function compute_rate!(dHdt, H, Hold, qH, dt, damp, dx, nx)
    Threads.@threads for ix=1:nx
        if (2<=ix<=nx-1)  dHdt[ix-1] = -(H[ix] - Hold[ix])/dt -(qH[ix]-qH[ix-1])/dx + damp*dHdt[ix-1]  end
    end
    return
end

function compute_update!(H, dHdt, dtau, nx)
    Threads.@threads for ix=1:nx
        if (2<=ix<=nx-1)  H[ix] = H[ix] + dtau*dHdt[ix-1]  end
    end
    return
end

@views function diffusion_1D()
    @show Threads.nthreads()
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
    qH   = zeros(nx-1)
    dHdt = zeros(nx-2)
    # Initial condition
    H    = exp.(.-(xc.-lx./2.0).^2)
    Hold = copy(H)
    H0   = copy(H)
    it = 1; err = 2*epsi
    # Time loop
    while err>epsi
        compute_flux!(qH, H, D, dx, nx)
        compute_rate!(dHdt, H, Hold, qH, dt, damp, dx, nx)
        compute_update!(H, dHdt, dtau, nx)
        it += 1; err = maximum(dHdt)
    end
    @printf("Total time = %1.2f, it tot = %d \n", round(dt, sigdigits=2), it)
    # Visualise
    plot(xc, H0, linewidth=3); display(plot!(xc, H, legend=false, framestyle=:box, linewidth=3, xlabel="lx", ylabel="H", title="damped diffusion (niter=$it)"))
    return
end

diffusion_1D()
