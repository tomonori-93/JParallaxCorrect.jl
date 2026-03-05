include("../src/JParallaxCorrect.jl")

using Test
using .JParallaxCorrect

d2r = π/180.0
r2d = 180.0/π

#-- Himawari-9 parameters
re = 6378.1370e3
rp = 6356.7523e3
hsat = 42164.0e3
lsat = 140.7
psat = 0.0
#-- Himawari-9 parameters

hmax = 15.0e3  # cloud top height
lonmin = 120.0
lonmax = 160.0
nlon = 4000
missing_value = -1.0
dlon = (lonmax - lonmin) / nlon
lon = lonmin .+ (Vector(1:nlon) .- 1) .* dlon
lon = reshape(lon,nlon,1)
lat = fill(0.0,nlon,1)
h_cld = fill(hmax,nlon,1)
lon2 = fill(0.0,nlon,2)
lat2 = fill(0.0,nlon,2)
h_cld2 = fill(hmax,nlon,2)
lon2[1:end,1] = lon[1:end,1]
lon2[1:end,2] = lon[1:end,1]
lat2[1:end,1] .= - dlon
lat2[1:end,2] .= dlon

@testset "parallax_correct_core" begin

    lon_grid, lat_grid = parallax_correct_core(lon.*d2r, lat.*d2r, h_cld,
                            re, rp, hsat, psat*d2r, lsat*d2r, missing_value)

    xs = (hsat / rp) * cos(psat*d2r) * cos(lsat*d2r)
    ys = (hsat / rp) * cos(psat*d2r) * sin(lsat*d2r)
    zs = (hsat / rp) * sin(psat*d2r)
    x0 = (re / rp) .* cos.(lat .* d2r) .* cos.(lon .* d2r)
    y0 = (re / rp) .* cos.(lat .* d2r) .* sin.(lon .* d2r)
    z0 = sin.(lat .* d2r)
    Xs = x0 .* zs .- z0 .* xs
    Ys = y0 .* zs .- z0 .* ys
    X0 = xs .- x0
    Y0 = ys .- y0
    Z0 = zs .- z0

    alpha = x0 .* X0 + y0 .* Y0
    beta = X0 .* X0 .+ Y0 .* Y0
    tpara = (alpha ./ beta) .* (-1.0 .+ sqrt.(1.0 .+ (beta ./ (alpha .* alpha)) .* (2.0 * re .* h_cld .+ h_cld .* h_cld) ./ (rp * rp)))
    xa = ((1.0 .- tpara) .* x0 .+ tpara .* xs) ./ (1.0 .+ h_cld ./ re)
    ya = ((1.0 .- tpara) .* y0 .+ tpara .* ys) ./ (1.0 .+ h_cld ./ re)

    lon_check = map((x,y) -> atan(y, x), xa, ya)
    @test lon_grid ≈ lon_check
end

@testset "ParallaxCorrect" begin

    println(lon2.*d2r,lat.*d2r)
    println(parallax_correct_core(lon2.*d2r, lat2.*d2r, h_cld2,
				  re, rp, hsat, psat*d2r, lsat*d2r, missing_value))
    h_grid, val_grid = ParallaxCorrect(lon2.*d2r, lat2.*d2r, h_cld2, h_cld2, 
                          lon[1:end,1].*d2r, lat[1:1,1].*d2r, re, rp, hsat, psat*d2r, lsat*d2r, missing_value)
    # val_grid: dummy (not used)

    xs = (hsat / rp) * cos(psat*d2r) * cos(lsat*d2r)
    ys = (hsat / rp) * cos(psat*d2r) * sin(lsat*d2r)
    zs = (hsat / rp) * sin(psat*d2r)
    x0 = (re / rp) .* cos.(lat .* d2r) .* cos.(lon .* d2r)
    y0 = (re / rp) .* cos.(lat .* d2r) .* sin.(lon .* d2r)
    z0 = sin.(lat .* d2r)
    Xs = x0 .* zs .- z0 .* xs
    Ys = y0 .* zs .- z0 .* ys
    X0 = xs .- x0
    Y0 = ys .- y0
    Z0 = zs .- z0

    alpha = x0 .* X0 + y0 .* Y0
    beta = X0 .* X0 .+ Y0 .* Y0
    tpara = (alpha ./ beta) .* (-1.0 .+ sqrt.(1.0 .+ (beta ./ (alpha .* alpha)) .* (2.0 * re .* h_cld .+ h_cld .* h_cld) ./ (rp * rp)))
    xa = ((1.0 .- tpara) .* x0 .+ tpara .* xs) ./ (1.0 .+ h_cld ./ re)
    ya = ((1.0 .- tpara) .* y0 .+ tpara .* ys) ./ (1.0 .+ h_cld ./ re)

    lon_check = map((x,y) -> atan(y, x), xa, ya)
    h_check = h_cld
    for i in 1:size(lon)[1]
        if lon_check[1,1] > lon[i,1] * d2r
            h_check[i,1] = missing_value
        end
	if lon_check[end,1] < lon[i,1] * d2r
            h_check[i,1] = missing_value
        end
    end

    println(h_grid[1:100,1])
    println(h_check[1:100,1])
    @test h_grid ≈ h_check
end

