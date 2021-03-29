
# julia-parallel-course-EGU21

#### [vEGU2021: SC4.6 Solving differential equations in parallel with Julia | Thu, 29 Apr, 16:00–17:00 (CEST)](https://meetingorganizer.copernicus.org/EGU21/session/38986)

👉 **Organisation notes:**
- 💡 All material presented during the short course will be uploaded here and made available to participants 5 days prior to the event on this repository.
- ⚠️ Live participation to this short course requires [registration to EGU's main meeting](https://egu21.eu/register.html).
- Due to the time-limited course schedule (60 min), an interactive overview will cover the course's [objectives](#objectives) replacing an extensive hands-on.
- Further interests in solving PDEs with Julia on GPUs❓
    - Sign-up (free) for a hands-on workshop at [JuliaCon2021](https://juliacon.org/2021/)
    - Check out this [online geo-HPC tutorial](https://github.com/luraess/geo-hpc-course)

----
This short course covers trendy areas in modern geocomputing with broad geoscientific applications. The physical processes governing natural systems' evolution are often mathematically described as systems of differential equations. A performant numerical implementation of the solving algorithm leveraging modern hardware is key and permits to tackle problems that were technically not possible a decade ago.


## Content
* [Objectives](#objectives)
* [Structure of the repository](#structure-of-the-repository)
* [Getting started](#getting-started) _(not part of the live course)_
* 👉 [Short course material](#short-course-material)
* [Extras](#extras) _(if time permits)_
* [Further reading](#further-reading)


## Objectives
The goal of this short course is to offer an interactive overview on how to solve systems of differential equations in parallel on GPUs using the [Julia language]. [Julia] combines high-level language simplicity and low-level language performance. The resulting codes and applications are fast, short and readable \[[1][JuliaCon20a], [2][JuliaCon20b], [3][JuliaCon19]\]. We will design and implement a numerical algorithm that predicts ice flow dynamics over mountainous topography (Greenland) using a high-performance computing approach:

![Greenland ice cap](docs/greenland_1.png)

The online course consists of 2 parts:
1. You will learn about the [Julia language], parallel and distributed computing and iterative solvers.
2. You will implement a PDE solver to predict ice flow dynamics on real topography.

By the end of this short course, you will:
- Have a GPU PDE solver that predicts ice-flow;
- Have a Julia code that achieves similar performance than legacy codes (C, CUDA, MPI);
- Know how the Julia language solves the "two-language problem";
- Be able to leverage the computing power of modern GPU accelerated servers and supercomputers.

⤴️ [_back to content_](#content)

## Structure of the repository
The course repository lists following folders and items:
- the [data](data) folder contains various low resolution Greenland input data (bedrock topography, surface elevation, ice thickness, masks, ...) downscaled from [BedMachine Greenland v3] - note the filenames include grid resolution information `(nx, ny)`;
- the [docs](docs) folder contains documentation linked in the [README](README.md);
- the output folder _will_ contain the various code output, mainly figures in png format;
- the [scripts](scripts) folder contains the scripts this course is about 🎉
- the [`Project.toml`](Project.toml) file is a Julia project file, tracking the used packages and enabling a reproducible environment.

⤴️ [_back to content_](#content)

## Getting started
> ⚠️ Due to the time limitation, the short course will not cover the [Getting started](#getting-started) steps. These are meant to provide directions to the participant willing to actively try out the examples during the short course. It is warmly recommended to perform the [Getting started](#getting-started) steps before the beginning of the workshop.

The provided directions will get you started with:
1. [Installing Julia v1.5.4](#installing-julia-v154) - 2 configurations are suggested:
- running Julia from the [terminal with an external text editor](#terminal--external-editor)
- running Julia from [VS Code](#vs-code)

2. [Running the scripts](#running-the-scripts) from the course repository.

> 👉 _**Note: This course relies on Julia v1.5.4**. The install configuration are tested on a MacBook Pro running macOS 10.15.7, a Linux GPU server running Ubuntu 20.04 LTS and a Linux GPU server running CentOS 8._

### Installing Julia v1.5.4
Check you have an active internet connexion and [download Julia v1.5.4](https://julialang.org/downloads/oldreleases/#v154_march_11_2021) for your platform following the install directions provided under [help] if needed.

Alternatively, open a terminal and download the binaries (select the one for your platform):
```sh
wget https://julialang-s3.julialang.org/bin/winnt/x64/1.5/julia-1.5.4-win64.exe # Windows
wget https://julialang-s3.julialang.org/bin/mac/x64/1.5/julia-1.5.4-mac64.dmg # macOS
wget https://julialang-s3.julialang.org/bin/linux/x64/1.5/julia-1.5.4-linux-x86_64.tar.gz # Linux x86
tar -xzf julia-1.5.4-<win64, mac64, linux-x86_64>.tar.gz # selecting your platform
```
Then add following lines in your `.bashrc`, `.profile`, or `config` file:
```sh
vim ~/.bashrc
PATH=<path-to>/julia-1.5.4/bin/:$PATH
export JULIA_CUDA_USE_BINARYBUILDER=false
```

#### Terminal + external editor
Ensure you have a text editor with syntax highlighting support for Julia. From within the terminal, type
```sh
julia
```
and you should be all set.

#### VS Code
If you'd enjoy a more IDE type of environment, [check out VS Code](https://code.visualstudio.com). Follow the [installation directions](https://github.com/julia-vscode/julia-vscode#getting-started) for the [Julia VS Code extension](https://www.julia-vscode.org).

### Running the scripts
To get started with the course,
1. clone (or download the ZIP archive) the course repository ([help here](https://docs.github.com/en/github/creating-cloning-and-archiving-repositories/cloning-a-repository))
```sh
git clone https://github.com/luraess/julia-parallel-course-EGU21.git
```
2. Navigate to `julia-parallel-course-EGU21`
```sh
cd julia-parallel-course-EGU21
```
3. From the terminal, launch Julia with the `--project` flag to read-in project environment related informations such as used packages
```sh
julia --project
```
3. From VS Code, follow the [instructions from the documentation](https://www.julia-vscode.org/docs/stable/gettingstarted/) to get started.

---
Now that you launched Julia, you should be in the [Julia REPL]. We now need to ensure all the packages we need to be installed before using them. To do so, enter the [Pkg mode](https://docs.julialang.org/en/v1/stdlib/REPL/#Pkg-mode) by typing `]`. Then, instantiate the project which should trigger the download of the packages. Exit the Pkg mode with CRTL+C:
```julia-repl
julia> ]

(julia-parallel-course-EGU21) pkg> st
Status `~/julia-parallel-course-EGU21/Project.toml`
  [4138dd39] JLD v0.12.2
  [85f8d34a] NCDatasets v0.11.3
  [94395366] ParallelStencil v0.5.1
  [91a5bcdd] Plots v1.11.1
  [37e2e46d] LinearAlgebra
  [de0858da] Printf

(julia-parallel-course-EGU21) pkg> instantiate
   Updating registry at `~/.julia/registries/General`
   Updating git-repo `https://github.com/JuliaRegistries/General.git`
   # [...]

julia>
```
To test your install, go to the [scripts](scripts) folder and run the [`iceflow.jl`](scripts/iceflow.jl) code. Note that you can execute shell commands from within the [Julia REPL] first typing `;`:
```julia-repl
julia> ;

shell> cd scripts/

julia> include("iceflow.jl")
```
Running this the first time will (pre-)complie the various installed packages and will take some time.  Subsequent runs, by executing `include("iceflow.jl")`, should take around 10s.

You should then see two figures saved in a newly created output folder, the second being the comparison between modelled and observed ice thickness distribution over Greenland:

![Greenland ice cap](docs/iceflow_out2.png)

### Multi-threading on CPUs
On the CPU, multi-threading is made accessible via [Base.Threads] and the environment variable [JULIA_NUM_THREADS] can be used to define the number of cores to use on the CPU, e.g. `export JULIA_NUM_THREADS=2` to enable 2 threads (2 CPU cores).

### Running on GPUs
The [CUDA.jl] module permits to launch compute kernels on Nvidia GPUs natively from within [Julia]. [JuliaGPU] provides further reading and [introductory material](https://juliagpu.gitlab.io/CUDA.jl/tutorials/introduction/) about GPU ecosystems within Julia.

⤴️ [_back to content_](#content)

## Short course material
This section lists the material discussed within this 60 min. short course:
* [Part 1 - Julia, parallel computing, iterative solvers](#part-1---julia-parallel-computing-iterative-solvers)
    * [Why Julia](#why-julia)
    * [Diffusion processes](#diffusion-processes)
    * [Iterative solvers](#iterative-solvers)
    * [Parallel GPU computing](#parallel-gpu-computing)
    * [XPU computing](#xpu-computing)
* [Part 2 - solving PDEs to predict ice flow](#part-2---solving-pdes-to-predict-ice-flow)
    * [SIA equation](#sia-equation)
    * [XPU SIA implementation](#xpu-sia-implementation)

💡 In this course we will implement a 2D nonlinear diffusion equation on GPUs in Julia using the finite-difference method and an iterative solving approach. The goal is to resolve the shallow ice approximation (SIA) and predict ice flow over Greenland.

⤴️ [_back to content_](#content)

---

### Part 1 - Julia, parallel computing, iterative solvers

#### Why Julia
_by M. Werder_

🚧 WIP

- the cool Julia ecosystem
- git integration
- modern code dev and CI
- ...

⤴️ [_back to course material_](#short-course-material)

#### Diffusion processes
Let's start with a simple 1D linear diffusion example to (i) see the difference between explicit and implicit and (ii) to compare the serial CPU vs the parallel GPU implementations. The diffusion of a quantity `H` over time `t` can be described as (1a) a diffusive flux, (1b) a flux balance and (1c) an update rule:
```md
qH    = -D*dH/dx  (1a)
dHdt  =  -dqH/dx  (1b)
dH/dt = dHdt      (1c)
```
The [`diffusion_1D_expl.jl`](scripts/diffusion_1D_expl.jl) code implements an iterative and explicit solution of eq. (1) for an initial Gaussian profile
```md
H = exp(-(x-lx/2.0)^2)
```

![](docs/diffusion_expl.png)

How to go with an implicit solution _**and**_ keeping it "matrix-free" ?

⤴️ [_back to course material_](#short-course-material)

#### Iterative solvers
The [`diffusion_1D_impl.jl`](scripts/diffusion_1D_impl.jl) code implements an iterative implicit solution of eq. (1). How ? We add the physical time derivative `dh/dt=(H-Hold)/dt` to the rate of change `dHdt` 
```md
dHdt = -(H-Hold)/dt -dqH/dx
```
and iterate until the values of `dHdt` (the residual of the eq. (1)) drop below a defined tolerance level `tol`.

![](docs/diffusion_impl.png)

It works, but the iteration count seems to be pretty high (`niter>1000`). There is a simple way to circumvent this by adding "damping" (`damp`) to the rate-of-change `dHdt`, analogous to adding friction to enable faster convergence
```md
dHdt = -(H-Hold)/dt -dqH/dx + damp*dHdt
```
The [`diffusion_1D_damp.jl`](scripts/diffusion_1D_damp.jl) code implements a damped iterative implicit solution of eq. (1). The iteration count drops to `niter<200`.

![](docs/diffusion_damp.png)

⤴️ [_back to course material_](#short-course-material)

#### Parallel GPU computing
So now we have a cool iterative and implicit solver in about 30 lines of code 🎉. Good enough for low resolution 1D calculations. What if we need more - 2D, 3D and high resolution to capture local and nonlinear physics ? Parallel and GPU computing makes it possible. Let's start from the [`diffusion_1D_damp.jl`](scripts/diffusion_1D_damp.jl) code and port it to GPU (with some intermediate steps).

1. Extract the physics calculations from [`diffusion_1D_damp.jl`](scripts/diffusion_1D_damp.jl), i.e. the time loop:
```julia
# [...] skipped lines
qH         .= -D*diff(H)/dx
dHdt       .= -(H[2:end-1].-Hold[2:end-1])/dt .-diff(qH)/dx .+ damp*dHdt
H[2:end-1] .= H[2:end-1] .+ dtau*dHdt
# [...] skipped lines
```

2. Split the calculations into separate functions (or kernels) and call those functions within the time loop. The [`diffusion_1D_damp_fun.jl`](scripts/diffusion_1D_damp_fun.jl) implements those modifications:
```julia
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
# [...] skipped lines
compute_flux!(qH, H, D, dx, nx)
compute_rate!(dHdt, H, Hold, qH, dt, damp, dx, nx)
compute_update!(H, dHdt, dtau, nx)
# [...] skipped lines
```
> 💡 Julia enables multi-threading capabilities by simply adding `Threads.@threads` to the outermost loop (here over `ix`).

3. The last step is to replace the (multi-threaded) loop by a vectorised index `ix = (blockIdx().x-1) * blockDim().x + threadIdx().x` specific to GPU execution. Each `ix` is executed concurrently on a different GPU thread. The [`diffusion_1D_damp_gpu.jl`](scripts/diffusion_1D_damp_gpu.jl) implements those modifications to run on GPUs:
```julia
using CUDA
# [...] skipped lines
function compute_flux!(qH, H, D, dx, nx)
    ix = (blockIdx().x-1) * blockDim().x + threadIdx().x
    if (ix<=nx-1)  qH[ix] = -D*(H[ix+1]-H[ix])/dx  end
    return
end

function compute_rate!(dHdt, H, Hold, qH, dt, damp, dx, nx)
    ix = (blockIdx().x-1) * blockDim().x + threadIdx().x
    if (2<=ix<=nx-1)  dHdt[ix-1] = -(H[ix] - Hold[ix])/dt -(qH[ix]-qH[ix-1])/dx + damp*dHdt[ix-1]  end
    return
end

function compute_update!(H, dHdt, dtau, nx)
    ix = (blockIdx().x-1) * blockDim().x + threadIdx().x
    if (2<=ix<=nx-1)  H[ix] = H[ix] + dtau*dHdt[ix-1]  end
    return
end
# [...] skipped lines
@cuda blocks=cublocks threads=cuthreads compute_flux!(qH, H, D, dx, nx)
synchronize()
@cuda blocks=cublocks threads=cuthreads compute_rate!(dHdt, H, Hold, qH, dt, damp, dx, nx)
synchronize()
@cuda blocks=cublocks threads=cuthreads compute_update!(H, dHdt, dtau, nx)
synchronize()
# [...] skipped lines
```
> 💡 We use `@cuda blocks=cublocks threads=cuthreads` to launch the GPU function on the appropriate number of threads, i.e. "parallel workers". The numerical grid resolution `nx` must now be chosen accordingly to the number of workers `nx=cublocks*cuthreads`.

⤴️ [_back to course material_](#short-course-material)

### XPU computing
Wouldn't it be great to be able to combine the multi-thread CPU and GPU implementations into a single "XPU" code to be able to run on various hardware with only changing a `USE_GPU` switch ? Using [ParallelStencil.jl] enables this, as well more other cool features. The [`diffusion_1D_damp_xpu.jl`](scripts/diffusion_1D_damp_xpu.jl) uses [ParallelStencil.jl] for an XPU implementation on various backends:

```julia
const USE_GPU = false
using ParallelStencil
using ParallelStencil.FiniteDifferences1D
@static if USE_GPU
    @init_parallel_stencil(CUDA, Float64, 1)
else
    @init_parallel_stencil(Threads, Float64, 1)
end
# [...] skipped lines
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
# [...] skipped lines
@parallel compute_flux!(qH, H, D, dx)
@parallel compute_rate!(dHdt, H, Hold, qH, dt, damp, dx)
@parallel compute_update!(H, dHdt, dtau)
# [...] skipped lines
```
Various macros `@(...)` permit to deal with the low-level technicalities and the `USE_GPU` flag enables to switch between CPU and GPU backend. The resulting code is short and readable and solves the "two-language problem"; development and production code implementations are regrouped into a single code. 

⤴️ [_back to course material_](#short-course-material)

### Part 2 - solving PDEs to predict ice flow

#### SIA equation
Let's move from the simple 1D linear diffusion example to the SIA equation, a 2D nonlinear diffusion equation:
```md
qHx   = -D*d(B+H)/dx              (2a)
qHy   = -D*d(B+H)/dy              (2b)
dHdt  = -(dqHx/dx + dqHy/dy) + M  (2c)
dH/dt = dHdt                      (2d)
```
where `B` is the bedrock elevation, `H` the ice thickness, `M` the mass balance (accumulation, ablation). The diffusion coefficient `D` is nonlinear and function of surface elevation `B+H` and the power-law exponent `n`:
```md
D = a*H^(npow+2)*sqrt((d(B+H)/dx)^2 + (d(B+H)/dy)^2)^(npow-1)
```
We implement climate forcing using a simple mass balance (accumulation, ablation) `M` formulation:
```md
M  = min(grad_b*(B+H - z_ELA), b_max)
```
as function of the surface elevation `B+H` and capped by the maximal accumulation rate `b_max`. The mass balance gradient `grad_b` is defined as
```md
grad_b = (1.3517 - 0.014158*LAT)/100.0*0.91
```
where `LAT` is the latitude (taken from \[[5][Machgut16]\]). The equilibrium line altitude (where accumulation = ablation) `z_ELA` is latitude dependent, ranging from 1300m (South) to 1000m (North) as suggsted by \[[5][Machgut16]\].

The [`iceflow.jl`](scripts/iceflow.jl) code implements the 2D SIA equations using the iterative implicit damped formulation as in [`diffusion_1D_damp.jl`](scripts/diffusion_1D_damp.jl). The calculation of the SIA PDEs resumes in these 13 lines of Julia code:
```julia
# mass balance
M     .= min.(grad_b.*(S .- z_ELA), b_max)
# compute diffusivity
dSdx  .= diff(S, dims=1)/dx
dSdy  .= diff(S, dims=2)/dy
gradS .= sqrt.(av_ya(dSdx).^2 .+ av_xa(dSdy).^2)
D     .= a*av(H).^(npow+2) .* gradS.^(npow-1)
# compute flux
qHx   .= .-av_ya(D).*diff(S[:,2:end-1], dims=1)/dx
qHy   .= .-av_xa(D).*diff(S[2:end-1,:], dims=2)/dy
# update ice thickness
dt    .= dtsc*min.(10.0, cfl./(epsi .+ av(D)))
ResH  .= .-(diff(qHx, dims=1)/dx .+ diff(qHy, dims=2)/dy) .+ inn(M)
dHdt  .= dHdt.*damp .+ ResH
H[2:end-1,2:end-1] .= max.(0.0, inn(H) .+ dt.*dHdt)
# apply mask
H[Mask.==0] .= 0.0
# update surface
S     .= B .+ H
```
🚧 WIP - to add:
- output figure
- some words on iteration count and time to solution on specific resolution (low res, slow)
- physical time derivative to be able to switch from steady state to specific time-span

#### XPU SIA implementation
Applying what we learned from the 1D diffusion equation, we can now instrument the [`iceflow.jl`](scripts/iceflow.jl) code (see code snippet just above) to make it XPU compatible using [ParallelStencil.jl]; the [`iceflow_xpu.jl`](scripts/iceflow_xpu.jl) code:
```julia
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
# [...] skipped lines
@parallel function compute_M_dS!(M, dSdx, dSdy, S, z_ELA, grad_b, b_max, dx, dy)
    @all(M)    = min(@all(grad_b)*(@all(S) - @all(z_ELA)), b_max)
    @all(dSdx) = @d_xa(S)/dx
    @all(dSdy) = @d_ya(S)/dy
    return
end

@parallel function compute_D!(D, H, dSdx, dSdy, a, npow)
    @all(D) = a*@pow(@av(H), (npow+2)) * @pow(sqrt(@av_ya(dSdx)*@av_ya(dSdx) + @av_xa(dSdy)*@av_xa(dSdy)), npow-1)
    return
end

@parallel function compute_qH_dt!(qHx, qHy, dt, D, S, dtsc, cfl, epsi, dx, dy)
    @all(qHx) = -@av_ya(D)*@d_xi(S)/dx
    @all(qHy) = -@av_xa(D)*@d_yi(S)/dy
    @all(dt)  = dtsc*min(10.0, cfl/(epsi + @av(D)))
    return
end

@parallel function compute_dHdt!(ResH, dHdt, qHx, qHy, M, damp, dx, dy)
    @all(ResH) = -(@d_xa(qHx)/dx + @d_ya(qHy)/dy) + @inn(M)
    @all(dHdt) = @all(dHdt)*damp + @all(ResH)
    return
end

@parallel function compute_H!(H, dHdt, dt)
    @inn(H) = max(0.0, @inn(H) + @all(dt)*@all(dHdt))
    return
end

@parallel_indices (ix,iy) function compute_Mask_S!(H, S, B, Mask)
    if (ix<=size(H,1) && iy<=size(H,2)) if (Mask[ix,iy]==0) H[ix,iy] = 0.0 end end
    if (ix<=size(H,1) && iy<=size(H,2)) S[ix,iy] = B[ix,iy] + H[ix,iy] end    
    return
end
# [...] skipped lines
@parallel compute_M_dS!(M, dSdx, dSdy, S, z_ELA, grad_b, b_max, dx, dy)
@parallel compute_D!(D, H, dSdx, dSdy, a, npow)
@parallel compute_qH_dt!(qHx, qHy, dt, D, S, dtsc, cfl, epsi, dx, dy)
@parallel compute_dHdt!(ResH, dHdt, qHx, qHy, M, damp, dx, dy)
@parallel compute_H!(H, dHdt, dt)
@parallel compute_Mask_S!(H, S, B, Mask)
```
🚧 Needs some words on high res example and maybe zoomed in figure into specific location.


⤴️ [_back to course material_](#short-course-material)

## Extras

🚧 WIP

### Performance metric
Curious about relevant performance metric for memory-bounded applications? Check out the [performance metric section](https://github.com/omlins/ParallelStencil.jl#performance-metric) from the [ParallelStencil.jl] module and this [JuliaCon2020][JuliaCon20a] presentation \[[1][JuliaCon20a]\].

Add a graph showing iteration count, time to solution, as fn of resolution ?...

### Multi-XPU implementation


### Simple inversion


⤴️ [_back to content_](#content)


## Further reading
\[1\] [Omlin, S., Räss, L., Kwasniewski, G., Malvoisin, B., & Podladchikov, Y. Y. (2020). Solving Nonlinear Multi-Physics on GPU Supercomputers with Julia. JuliaCon Conference, virtual.][JuliaCon20a]

\[2\] [Räss, L., Reuber, G., Omlin, S. (2020). Multi-Physics 3-D Inversion on GPU Supercomputers with Julia. JuliaCon Conference, virtual.][JuliaCon20b]

\[3\] [Räss, L., Omlin, S., & Podladchikov, Y. Y. (2019). Porting a Massively Parallel Multi-GPU Application to Julia: a 3-D Nonlinear Multi-Physics Flow Solver. JuliaCon Conference, Baltimore, USA.][JuliaCon19]

\[4\] [Frankel, S. P. (1950). Convergence rates of iterative treatments of partial differential equations, Mathe. Tables Other Aids Comput., 4, 65–75.][Frankel50]

\[5\] [Machgut, H. et al. (2016). Greenland surface mass-balance observations from the ice-sheet ablation area and local glaciers. Journal of Glaciology, 62(235), 861-887.][Machgut16]

⤴️ [_back to content_](#content)


[Julia]: https://julialang.org
[Julia language]: https://docs.julialang.org/en/v1/
[Julia REPL]: https://docs.julialang.org/en/v1/stdlib/REPL/
[Base.Threads]: https://docs.julialang.org/en/v1/base/multi-threading/
[JULIA_NUM_THREADS]:https://docs.julialang.org/en/v1.0.0/manual/environment-variables/#JULIA_NUM_THREADS-1
[CUDA.jl]: https://github.com/JuliaGPU/CUDA.jl
[JuliaGPU]: https://juliagpu.org
[ParallelStencil.jl]: https://github.com/omlins/ParallelStencil.jl

[BedMachine Greenland v3]: https://sites.uci.edu/morlighem/dataproducts/bedmachine-greenland/

[JuliaCon20a]: https://www.youtube.com/watch?v=vPsfZUqI4_0
[JuliaCon20b]: https://www.youtube.com/watch?v=1t1AKnnGRqA
[JuliaCon19]: https://www.youtube.com/watch?v=b90qqbYJ58Q
[Frankel50]: /docs/frankel_1950.pdf
[Machgut16]: https://doi.org/10.1017/jog.2016.75
