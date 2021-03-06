# Solver better capabale of running high-resolution simulations, using a few more tricks.
const USE_GPU = true
using ParallelStencil
using ParallelStencil.FiniteDifferences2D
@static if USE_GPU
    @init_parallel_stencil(CUDA, Float64, 2)
else
    @init_parallel_stencil(Threads, Float64, 2)
end
using JLD, Plots, Printf, LinearAlgebra

# CPU functions
@views av(A)  = 0.25*(A[1:end-1,1:end-1].+A[2:end,1:end-1].+A[1:end-1,2:end].+A[2:end,2:end])
@views inn(A) = A[2:end-1,2:end-1]

# GPU function
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

@parallel function compute_M_dS!(M, dSdx, dSdy, S, z_ELA, grad_b, b_max, dx, dy)
    @all(M)    = min(@all(grad_b)*(@all(S) - @all(z_ELA)), b_max)
    @all(dSdx) = @d_xa(S)/dx
    @all(dSdy) = @d_ya(S)/dy
    return
end

@parallel function compute_D!(D, H, dSdx, dSdy, a, npow)
    @all(D) = a*@av(H)^(npow+2) * sqrt(@av_ya(dSdx)*@av_ya(dSdx) + @av_xa(dSdy)*@av_xa(dSdy))^(npow-1)
    return
end

@parallel function compute_flux!(qHx, qHy, D, S, dx, dy)
    @all(qHx)  = -@av_ya(D)*@d_xi(S)/dx
    @all(qHy)  = -@av_xa(D)*@d_yi(S)/dy
    return
end

@parallel function compute_dHdt!(dtau, ResH, dHdt, D, qHx, qHy, M, dtausc, cfl, epsi, damp, dx, dy)
    @all(dtau) = dtausc*min(10.0, cfl/(epsi + @av(D)))
    @all(ResH) = -(@d_xa(qHx)/dx + @d_ya(qHy)/dy) + @inn(M)
    @all(dHdt) = @all(dHdt)*damp + @all(ResH)
    return
end

@parallel function compute_minloc!(A2, A)
    @inn(A2) = @minloc(A)
    return
end

@parallel function compute_H!(H, dHdt, dtau)
    @inn(H) = max(0.0, @inn(H) + @all(dtau)*@all(dHdt))
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

@views function iceflow(dx, dy, Zbed, Hice, Mask, grad_b, z_ELA, b_max)
    println("Initialising ice flow model ... ")
    # physics
    s2y      = 3600*24*365.25  # seconds to years
    rho_i    = 910.0           # ice density
    g        = 9.81            # gravity acceleration
    npow     = 3.0             # Glen's power law exponent
    a0       = 1.5e-24         # Glen's law enhancement term
    # numerics
    @assert (dx>0 && dy>0) "dx and dy need to be positive"
    nx, ny   = size(Zbed,1), size(Zbed,2) # numerical grid resolution
    @assert (nx, ny) == size(Zbed) == size(Hice) == size(Mask) "Sizes don't match"
    itMax    = 1e5             # number of iteration (max)
    nout     = 1000            # error check frequency
    tolnl    = 5e-4            # nonlinear tolerance
    epsi     = 1e-4            # small number
    damp     = 0.72             # convergence accelerator
    dtausc   = 1.0/4.0         # iterative dtau scaling
    ntloc    = 2               # local min iterations
    # derived physics
    a        = 2.0*a0/(npow+2)*(rho_i*g)^npow*s2y
    # derived numerics
    cfl      = max(dx^2,dy^2)/4.1
    # array initialization
    Err      = @zeros(nx  , ny  )
    dSdx     = @zeros(nx-1, ny  )
    dSdy     = @zeros(nx  , ny-1)
    D        = @zeros(nx-1, ny-1)
    qHx      = @zeros(nx-1, ny-2)
    qHy      = @zeros(nx-2, ny-1)
    dtau     = @zeros(nx-2, ny-2)
    dtau2    = @zeros(nx-2, ny-2)
    ResH     = @zeros(nx-2, ny-2)
    dHdt     = @zeros(nx-2, ny-2)
    Vx       = @zeros(nx-1, ny-1)
    Vy       = @zeros(nx-1, ny-1)
    M        = @zeros(nx  , ny  )
    S        = @zeros(nx  , ny  )
    grad_b   = Data.Array(grad_b)
    z_ELA    = Data.Array(z_ELA)
    Mask     = Data.Array(Mask)
    # initial conditions
    B        = Data.Array(Zbed)
    H        = Data.Array(Hice)
    S       .= B .+ H
    println(" starting iteration loop:")
    # iteration loop
    iter = 1; err = 2*tolnl
    while err>tolnl && iter<itMax
        @parallel compute_Err1!(Err, H)
        @parallel compute_M_dS!(M, dSdx, dSdy, S, z_ELA, grad_b, b_max, dx, dy)
        @parallel compute_D!(D, H, dSdx, dSdy, a, npow)
        @parallel compute_flux!(qHx, qHy, D, S, dx, dy)
        @parallel compute_dHdt!(dtau, ResH, dHdt, D, qHx, qHy, M, dtausc, cfl, epsi, damp, dx, dy)
        for iloc = 1:ntloc
            @parallel compute_minloc!(dtau2, dtau)
            @parallel (1:size(dtau2,2)) bc_x!(dtau2)
            @parallel (1:size(dtau2,1)) bc_y!(dtau2)
            dtau, dtau2 = dtau2, dtau
        end
        @parallel compute_H!(H, dHdt, dtau)
        @parallel compute_Mask_S!(H, S, B, Mask)
        # error check
        if mod(iter, nout)==0
            @parallel compute_Err2!(Err, H)
            err = norm(Err)/length(Err)
            @printf(" iter = %d, error = %1.2e \n", iter, err)
            if isnan(err)
                error("""NaNs encountered.  Try a combination of:
                             decreasing `damp` and/or `dtausc`, more smoothing steps""")
            end
        end
        iter += 1
    end
    @parallel compute_Vel!(Vx, Vy, D, H, dSdx, dSdy, epsi)
    # return as GeoArrays
    return  as_geoarray(Array(H),  Zbed, name=:thickness),
            as_geoarray(Array(S),  Zbed, name=:surface),
            as_geoarray(Array(M),  Zbed, name=:smb),
            as_geoarray(Array(Vx), Zbed, name=:vel_x, staggerd=true),
            as_geoarray(Array(Vy), Zbed, name=:vel_y, staggerd=true)
end
# ------------------------------------------------------------------------------
include("../scripts/helpers.jl")

# load the data
print("Loading the data ... ")
Zbed, Hice, Mask, dx, dy, xc, yc = load_data(; nx=96) # nx=96,160 are included in the repo
                                                                      # other numbers will trigger a 2GB download
println("done.")

# apply some smoothing
print("Applying some smoothing ... ")
for is=1:2 # two smoothing steps
    smooth!(Zbed)
    smooth!(Hice)
end
println("done.")

# calculate mass balance coefficients for given spatial grid
grad_b, z_ELA, b_max = mass_balance_constants(xc, yc)

# run the SIA flow model
H, S, M, Vx, Vy = iceflow(dx, dy, Zbed, Hice, Mask, grad_b, z_ELA, b_max)

# handle output
do_visu = true
do_save = true

# visu and save
nx, ny = size(H)
if do_visu
    !ispath("../output") && mkdir("../output")

    H_v = copy(H); H_v[Mask.==0].=NaN
    Hice_v = copy(Hice); Hice_v[Mask.==0].=NaN
    S_v = copy(S); S_v[Mask.==0].=NaN
    M_v = copy(M); M_v[Mask.==0].=NaN
    V_v = sqrt.(Vx.^2 .+ Vy.^2)

    # outputs
    fontsize  = 7
    opts = (aspect_ratio=1, yaxis=font(fontsize, "Courier"), xaxis=font(fontsize, "Courier"),
            ticks=nothing, framestyle=:box, titlefontsize=fontsize, titlefont="Courier", colorbar_title="",
            xlabel="", ylabel="", xlims=(dims(H_v)[1][1],dims(H_v)[1][end]), ylims=(dims(H_v)[2][end],dims(H_v)[2][1]),
            )
    p1 = heatmap(S_v; c=:davos, title="Surface elev. [m]", opts...)
    p2 = heatmap(H_v; c=:davos, title="Ice thickness [m]", opts...)
    p3 = heatmap(log10.(V_v); clims=(0.1, 2.0), title="log10(vel) [m/yr]", opts...)
    p4 = heatmap(M_v; c=:devon, title="Mass Bal. rate [m/yr]", opts...)
    p = plot(p1, p2, p3, p4, size=(400,400), dpi=200) #background_color=:transparent, foreground_color=:white
    ## uncomment if you want a pop-up plot pane showing:
    # display(p)
    savefig("../output/iceflow_out1_xpu_$(nx)x$(ny).png")

    # error
    H_diff = Hice_v.-H_v
    fontsize = 7
    p1 = heatmap(Hice_v; c=:davos, title="H data [m]", opts...)
    p2 = heatmap(H_v; c=:davos, title="H model [m]", opts...)
    clim = max(abs.(extrema(H_diff[.!isnan.(H_diff)]))...)
    p3 = heatmap(H_diff; title="H (data-model) [m]",  clims=(-clim,clim), seriescolor=:balance, opts...)
    p = plot(p1, p2, p3, layout=(1, 3), size=(500,160), dpi=200) #background_color=:transparent, foreground_color=:white
    ## uncomment if you want a pop-up plot pane showing:
    # display(p)
    savefig("../output/iceflow_out2_xpu_$(nx)x$(ny).png")
end

if do_save
    save("../output/iceflow_xpu_$(nx)x$(ny).jld", "Hice", Hice,
                                                  "Mask", Mask,
                                                  "H"   , H,
                                                  "S"   , S,
                                                  "M"   , M,
                                                  "Vx"  , Vx,
                                                  "Vy"  , Vy,
                                                  "xc", xc, "yc", yc)
end

println("... done.")
