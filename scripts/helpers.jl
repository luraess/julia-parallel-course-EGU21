using Downloads, GeoData, NCDatasets, JLD

"""
    smooth!(A)

Smooth data contained in a matrix with one time step (CFL) of diffusion.
"""
@views function smooth!(A)
    A[2:end-1,2:end-1] .= A[2:end-1,2:end-1] .+ 1.0./4.1.*(diff(diff(A[:,2:end-1], dims=1), dims=1) .+ diff(diff(A[2:end-1,:], dims=2), dims=2))
    A[1,:]=A[2,:]; A[end,:]=A[end-1,:]; A[:,1]=A[:,2]; A[:,end]=A[:,end-1]
    return
end

datadir = "../data"
bm_file = joinpath(datadir, "BedMachineGreenland-2017-09-20.nc")

"""
    download_bedmachine_greenland(; force=false, force_earthdata=false, verbose=true)

Try to download the Bedmachine Greenland dataset.  During the short-course
from a server of ours.  Outside the course, it is behind a password-wall and
it is a bit more complicated, try and see the error message.

It will work if:
- there are credentials in the `~/.netrc` file for `urs.earthdata.nasa.gov`, see [1]
- or *during the course* via the copy we'll provide not needing a login

Keyword-arguments:
- set `force==true` to download even if the file is present.

Returns the local path to the file.

Refs:
[1] https://nsidc.org/support/how/how-do-i-programmatically-request-data-services
[2] https://github.com/JuliaLang/Downloads.jl/issues/103#issuecomment-818939509
"""
function download_bedmachine_greenland(; force=false, force_earthdata=false, verbose=true)
    verbose && @warn "This will download the 2GB file BedMachineGreenland-2017-09-20.nc.  Hit Ctrl-C to abort."
    url_course = "https://people.ee.ethz.ch/~werderm/6458-gpu-course/BedMachineGreenland-2017-09-20.nc"
    url = "https://n5eil01u.ecs.nsidc.org/ICEBRIDGE/IDBMG4.003/1993.01.01/BedMachineGreenland-2017-09-20.nc"

    if isfile(bm_file) && filesize(bm_file)==2249644927 && !force
        @warn "File exists and is of right size, not downloading."
        return bm_file
    end

    # try the uploaded file on my web-server
    if !force_earthdata
        try
            println("") # the _download_progress deletes a line, so provide one
            Downloads.download(url_course, bm_file; progress=_download_progress)
            return bm_file
        catch e
            if !(e isa Downloads.RequestError)
                println("Encountered an unexpected error but will try different source.  Error:\n $e")
            end
        end
    end

    # try earthdata/NASA
    # uses https://github.com/JuliaLang/Downloads.jl/issues/103#issuecomment-818939509
    # which may stop working at some point...
    downloader = Downloads.Downloader()
    easy_hook = (easy, info) -> begin
        Downloads.Curl.setopt(easy, Downloads.Curl.CURLOPT_NETRC, Downloads.Curl.CURL_NETRC_OPTIONAL)
        Downloads.Curl.setopt(easy, Downloads.Curl.CURLOPT_COOKIEFILE, "")
    end
    downloader.easy_hook = easy_hook

    try
        Downloads.download(url, out; downloader=downloader, progress=_download_progress)
        return bm_file
    catch e
        if !(e isa Downloads.RequestError)
            println("Encountered an unexpected error:\n $e")
        end
    end

    # no luck, error:
    error("""
        Could not download the Greenland Bedmachine data.  This can have a few reasons:
        - you don't have (correct) credentials to https://urs.earthdata.nasa.gov/home in
          the file `~/.netrc`.
        - something else...

        Either way, it is probably easiest if you download the file by hand (you will still
        need to create a `earthdata.nasa.gov` login).  To do this, go to
        https://nsidc.org/data/idbmg4 and click "Login to Earthdata".  Once logged in visit
        https://n5eil01u.ecs.nsidc.org/ICEBRIDGE/IDBMG4.003/1993.01.01/BedMachineGreenland-2017-09-20.nc
        and save the file.

        Move the file `BedMachineGreenland-2017-09-20.nc` to the folder
        $(joinpath(@__DIR__, "../data")).
        """)
end

"""
    load_bedmachine_greenland(;downscale=nothing, nx=nothing)

Load the Bedmachine data into memory at a specific resolution.  The resolution can be
chosen by setting either:
- the downscaling factor, between 1 and 329, or
- nx, between 10208 and 32 [default=96]

Note that this will automatically produce a dataset with x-grid
divisible by 32 and y-grid by 8 (as suitable for the GPU), and thus the exact
dimensions/downscaling maybe different.

Return:
- Zbed, Hice, Mask, dx, dy
"""
function load_bedmachine_greenland(;downscale=nothing, nx=96)
    if downscale!=nothing && nx!=96
        error("Only choose of the input augments `downscale` and `nx`")
    end
    # GPU is most performant for these divisors of the array-size
    nx_div = 32
    ny_div = 8

    res = get_resolution(;ds=downscale, nx=nx)
    nx = res.nx32

    if nx in [96, 160] # get from jld files in the repo
        data = if nx == 96
            load("../data/BedMachineGreenland_96_184.jld") # ultra low res data
        else
            load("../data/BedMachineGreenland_160_304.jld") # low res data
        end
        Hice, Mask, Zbed = data["Hice"], data["Mask"], data["Zbed"]
        xc, yc, dx, dy   = data["xc"], data["yc"], data["dx"], data["dy"]
    else  # get from nc file
        if !(isfile(bm_file) && filesize(bm_file)==2249644927)
            println("""Downloading the Bedmachine Greenland dataset.  This may take a while (2GB in size).
                          Hit Ctrl+C to abort.
                    """)
            download_bedmachine_greenland(verbose=false)
        end

        bms = NCDstack(bm_file)
        Zbed = downscale_and_crop(bms[:bed], res)
        Hice = downscale_and_crop(bms[:thickness], res)
        Mask = downscale_and_crop(bms[:mask], res)
        # mask: 0=ocean, 1=ice-free land, 2=grounded-ice, 3=floating-ice, 4=non-Greenland land
        Mask = (Mask.==1) .| (Mask.==2)
        dx, dy = abs.(step.(dims(Zbed)))
    end

    return Zbed, Hice, Mask, dx, dy
end


"Enumeration of ideal resolutions for different downscalings, nx, and ny"
resolutions =
[(ds = 1, nx32 = 10208, cropx = 10, ny8 = 18344, cropy = 2)
 (ds = 2, nx32 = 5088, cropx = 43, ny8 = 9168, cropy = 11)
 (ds = 3, nx32 = 3392, cropx = 44, ny8 = 6112, cropy = 12)
 (ds = 4, nx32 = 2528, cropx = 109, ny8 = 4584, cropy = 13)
 (ds = 5, nx32 = 2016, cropx = 142, ny8 = 3664, cropy = 30)
 (ds = 6, nx32 = 1696, cropx = 47, ny8 = 3056, cropy = 15)
 (ds = 7, nx32 = 1440, cropx = 144, ny8 = 2616, cropy = 40)
 (ds = 8, nx32 = 1248, cropx = 241, ny8 = 2288, cropy = 49)
 (ds = 9, nx32 = 1120, cropx = 146, ny8 = 2032, cropy = 66)
 (ds = 10, nx32 = 992, cropx = 307, ny8 = 1832, cropy = 35)
 (ds = 11, nx32 = 928, cropx = 20, ny8 = 1664, cropy = 52)
 (ds = 12, nx32 = 832, cropx = 245, ny8 = 1528, cropy = 21)
 (ds = 13, nx32 = 768, cropx = 246, ny8 = 1408, cropy = 54)
 (ds = 14, nx32 = 704, cropx = 375, ny8 = 1304, cropy = 103)
 (ds = 15, nx32 = 672, cropx = 152, ny8 = 1224, cropy = 0)
 (ds = 16, nx32 = 608, cropx = 505, ny8 = 1144, cropy = 57)
 (ds = 17, nx32 = 576, cropx = 442, ny8 = 1080, cropy = 2)
 (ds = 18, nx32 = 544, cropx = 443, ny8 = 1016, cropy = 75)
 (ds = 19, nx32 = 512, cropx = 508, ny8 = 960, cropy = 124)
 (ds = 21, nx32 = 480, cropx = 158, ny8 = 872, cropy = 54)
 (ds = 22, nx32 = 448, cropx = 383, ny8 = 832, cropy = 63)
 (ds = 24, nx32 = 416, cropx = 257, ny8 = 760, cropy = 129)
 (ds = 26, nx32 = 384, cropx = 259, ny8 = 704, cropy = 67)
 (ds = 29, nx32 = 352, cropx = 38, ny8 = 632, cropy = 46)
 (ds = 32, nx32 = 320, cropx = 9, ny8 = 568, cropy = 201)
 (ds = 35, nx32 = 288, cropx = 172, ny8 = 520, cropy = 180)
 (ds = 40, nx32 = 256, cropx = 17, ny8 = 456, cropy = 145)
 (ds = 45, nx32 = 224, cropx = 182, ny8 = 408, cropy = 30)
 (ds = 53, nx32 = 192, cropx = 94, ny8 = 344, cropy = 166)
 (ds = 63, nx32 = 160, cropx = 200, ny8 = 288, cropy = 264)
 (ds = 79, nx32 = 128, cropx = 184, ny8 = 232, cropy = 96)
 (ds = 104, nx32 = 96, cropx = 337, ny8 = 176, cropy = 145)
 (ds = 162, nx32 = 64, cropx = 11, ny8 = 112, cropy = 363)
 (ds = 329, nx32 = 32, cropx = 18, ny8 = 56, cropy = 250)]

"""
     get_resolution(;ds=nothing, nx=nothing, ny=nothing)

Return the resolution for the BedMachine Greenland dataset which
fits best with one of the desired
- downscaling (ds),
- number of x-gridpoints (nx), or
- number of y-gridpoints (ny)
such that the grid is of a good shape for performant GPU code.
Namely nx is divisible by 32 and ny by 8.

Downscale and crop the BedMachine data like so

    downscale_and_crop(A, get_resolution(...))
"""
function get_resolution(;ds=nothing, nx=nothing, ny=nothing)
    @assert sum([ds,nx,ny].==nothing)==2 "Only one option ds, nx, or ny can be used at once."
    i = if ds!=nothing
        findclosest(r -> abs(r.ds-ds), resolutions)
    elseif nx!=nothing
        findclosest(r -> abs(r.nx32-nx), resolutions)
    elseif ny!=nothing
        findclosest(r -> abs(r.ny8-ny), resolutions)
    end
    return resolutions[i]
end

"""
    findclosest(fn, collection)

Find item in `collection` for which `fn` returns
the lowest value. If there is a tie, pick the first item encountered.
Returns index.

    findclosest(x -> abs(3 - x), [0, 1, 6, 5]) # -> 2
"""
function findclosest(fn, collection)
    i = 0
    dist = typemax(fn(collection[1]))
    for (ii,it) in enumerate(collection)
        if fn(it) < dist
            dist = fn(it)
            i = ii
        end
    end
    return i
end

"Apply downscaling and cropping returned from `get_resolution`."
downscale_and_crop(A, res) = A[1 + floor(Int, res.cropx/2):res.ds:end - ceil(Int, res.cropx/2),
                               1 + floor(Int, res.cropy/2):res.ds:end - ceil(Int, res.cropy/2)]
