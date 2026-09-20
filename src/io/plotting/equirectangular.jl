#---------------------------------------------------------------------------------
# equirectangular.jl — spherical (λ, φ) fields onto a plate-carrée raster.
#
# Everything this file draws lives on the SHELL, i.e. on a 2-D manifold whose
# nodes carry (lon, lat) and nothing else. A contour plot of such a field needs
# a flat canvas, and the standard one — the one every shallow-water paper uses
# for the Galewsky jet, Williamson's test cases and the POD modes computed from
# them — is the EQUIRECTANGULAR (plate carrée) projection:
#
#     x = λ  [deg, -180 … 180] ,    y = φ  [deg, -90 … 90]
#
# i.e. longitude and latitude used directly as Cartesian plot coordinates. It is
# neither conformal nor equal-area; it is used because it is the identity map on
# the coordinates the data already carries, so nothing in the picture is an
# artefact of the projection.
#
# HOW THE RASTER IS BUILT, and why not by interpolation of scattered points.
# The nodes are the LGL nodes of a spectral element grid: they are NOT scattered
# — they come with the element connectivity, which tiles the sphere exactly.
# Rasterizing that tiling (each (ngl-1)² sub-quad of each element split into two
# triangles, filled by barycentric interpolation) is therefore a RENDERING of
# the discrete solution, not a re-interpolation of it:
#
#   * no smoothing parameter and no neighbour count to tune, hence nothing that
#     can invent or erase a feature — the failure mode of inverse-distance or
#     spline fits of the same data, which round off exactly the small-scale
#     structure the higher POD modes consist of;
#   * values are convex combinations of nodal values, so the raster cannot
#     overshoot the data. A mode plotted this way has the same extrema as the
#     mode itself, which is what makes the colour scale meaningful;
#   * the panel seams of the cubed sphere are invisible, because the triangles
#     are drawn from the connectivity and meet exactly there.
#
# THE TWO PLACES THE PROJECTION IS SINGULAR are handled explicitly:
#
#   * the DATELINE. λ = atan(y,x) jumps by 2π across it, so a sub-quad straddling
#     it would be drawn as a band stretching the whole way round the map. Each
#     quad's corner longitudes are therefore UNWRAPPED onto a branch containing
#     its first corner, and the quad is then drawn at λ, λ-360 and λ+360; only
#     the copies that overlap the canvas cost anything. A feature crossing the
#     dateline comes out continuous, entering at +180 and leaving at -180.
#
#   * the POLES, where the map is genuinely degenerate: a whole neighbourhood of
#     the pole collapses onto the single line φ = ±90, and the longitude of a
#     node sitting exactly on the pole (the cubed sphere puts one there whenever
#     the panel has an even number of elements per edge) is undefined —
#     atan(0,0) returns 0. Triangles touching the pole are slivers in λ and
#     leave a few pixels of the top and bottom rows uncovered. Those are filled
#     from the NEAREST NODE IN 3-D (chordal distance on the sphere, where there
#     is no singularity), which is the one choice that cannot depend on the
#     arbitrary λ of the polar node.
#
# S. Marras & contributors
#---------------------------------------------------------------------------------

export equirectangular_grid, equirectangular_raster


"""
    equirectangular_grid(nlon, nlat) -> (λ, φ)

Pixel CENTRES of an `nlon × nlat` plate-carrée canvas, in degrees:

    λ[i] = -180 + (i-½)·360/nlon ,   φ[j] = -90 + (j-½)·180/nlat

Centres rather than edges because `Plots.contourf(λ, φ, F')` reads its first two
arguments as the coordinates AT WHICH `F` is sampled.
"""
function equirectangular_grid(nlon::Int, nlat::Int)
    nlon > 1 && nlat > 1 ||
        error(" # ERROR equirectangular.jl: need nlon > 1 and nlat > 1, got ($nlon, $nlat).")
    dλ = 360.0/nlon
    dφ = 180.0/nlat
    λ  = [-180.0 + (i - 0.5)*dλ for i = 1:nlon]
    φ  = [ -90.0 + (j - 0.5)*dφ for j = 1:nlat]
    return λ, φ
end


"""
    equirectangular_raster(f, lon, lat, connijk, nelem, ngl; kwargs...) -> (λ, φ, F)

Render the nodal field `f` (length `npoin`) onto an equirectangular raster.
`lon`/`lat` are the nodal spherical coordinates IN RADIANS (`mesh.lon`,
`mesh.lat`), and `connijk[iel,i,j]` the spectral element connectivity.

Returns the pixel-centre coordinate vectors (degrees) and the raster `F`, of
size `nlon × nlat` — index order `(λ, φ)`, so pass `F'` to `Plots.contourf`,
as the flat-case plotters in jeplots.jl do with their own rasters.

Keyword arguments:

  * `nlon`, `nlat`  — raster size; the default 720 × 360 is ½° and oversamples
                      the shipped cubed sphere (~1.4° between nodes) about
                      threefold, which is what keeps the triangle edges from
                      showing.
  * `lfill_gaps`    — fill pixels no triangle covered (the polar caps, see the
                      header) from the nearest node. `false` leaves them `NaN`,
                      which `Plots` draws as blank.
  * `nfill_sample`  — how many nodes the nearest-node fallback searches. The
                      full node set is subsampled to roughly this many, since
                      the fallback is brute force and only runs on the handful
                      of uncovered pixels.
"""
function equirectangular_raster(f::AbstractVector, lon::AbstractVector, lat::AbstractVector,
                                connijk, nelem::Int, ngl::Int;
                                nlon::Int = 720, nlat::Int = 360,
                                lfill_gaps::Bool = true,
                                nfill_sample::Int = 4000)

    npoin = min(length(f), length(lon), length(lat))
    npoin > 0 || error(" # ERROR equirectangular.jl: empty field handed to the rasterizer.")

    λ, φ = equirectangular_grid(nlon, nlat)
    F    = fill(NaN, nlon, nlat)

    dλ = 360.0/nlon
    dφ = 180.0/nlat
    λ0 = -180.0
    φ0 =  -90.0

    rad2deg_ = 180.0/π

    # corner buffers, reused by every sub-quad
    xq = zeros(Float64, 4)
    yq = zeros(Float64, 4)
    vq = zeros(Float64, 4)

    @inbounds for iel = 1:nelem
        for j = 1:ngl-1, i = 1:ngl-1

            ip1 = connijk[iel, i,   j  ]
            ip2 = connijk[iel, i+1, j  ]
            ip3 = connijk[iel, i+1, j+1]
            ip4 = connijk[iel, i,   j+1]

            (ip1 in 1:npoin && ip2 in 1:npoin && ip3 in 1:npoin && ip4 in 1:npoin) || continue

            vq[1] = Float64(f[ip1]); vq[2] = Float64(f[ip2])
            vq[3] = Float64(f[ip3]); vq[4] = Float64(f[ip4])
            (isfinite(vq[1]) && isfinite(vq[2]) && isfinite(vq[3]) && isfinite(vq[4])) || continue

            yq[1] = Float64(lat[ip1])*rad2deg_; yq[2] = Float64(lat[ip2])*rad2deg_
            yq[3] = Float64(lat[ip3])*rad2deg_; yq[4] = Float64(lat[ip4])*rad2deg_

            #
            # UNWRAP onto the branch of corner 1. Without this a quad with
            # corners at λ = 179.5 and λ = -179.5 is drawn as a 359°-wide band.
            #
            xq[1] = Float64(lon[ip1])*rad2deg_
            for c = 2:4
                ipc  = c == 2 ? ip2 : (c == 3 ? ip3 : ip4)
                xc   = Float64(lon[ipc])*rad2deg_
                while xc - xq[1] >  180.0; xc -= 360.0; end
                while xc - xq[1] < -180.0; xc += 360.0; end
                xq[c] = xc
            end

            xmin = min(xq[1], xq[2], xq[3], xq[4])
            xmax = max(xq[1], xq[2], xq[3], xq[4])

            #
            # …and draw the unwrapped quad at every 360° offset that still
            # touches the canvas, so that what leaves at +180 re-enters at -180.
            #
            for shift in (-360.0, 0.0, 360.0)
                (xmax + shift < λ0 || xmin + shift > λ0 + 360.0) && continue
                # two triangles, counter-clockwise as the sub-quad was built
                _raster_triangle!(F, xq[1]+shift, yq[1], vq[1],
                                     xq[2]+shift, yq[2], vq[2],
                                     xq[3]+shift, yq[3], vq[3],
                                  λ0, φ0, dλ, dφ, nlon, nlat)
                _raster_triangle!(F, xq[1]+shift, yq[1], vq[1],
                                     xq[3]+shift, yq[3], vq[3],
                                     xq[4]+shift, yq[4], vq[4],
                                  λ0, φ0, dλ, dφ, nlon, nlat)
            end
        end
    end

    lfill_gaps && _fill_raster_gaps!(F, λ, φ, f, lon, lat, npoin, nfill_sample)

    return λ, φ, F
end


#
# Convenience overload: everything the rasterizer needs is on the mesh.
# Deliberately UNTYPED in `mesh` so that this file stays free of Jexpresso types
# and can be included on its own (test/test_pod.jl does exactly that).
#
function equirectangular_raster(f::AbstractVector, mesh; kwargs...)
    return equirectangular_raster(f, mesh.lon, mesh.lat, mesh.connijk,
                                  Int(mesh.nelem), Int(mesh.ngl); kwargs...)
end


#
# Fill one triangle by barycentric interpolation. Pixel centres only: a pixel is
# painted when its CENTRE falls inside, which is what makes two triangles
# sharing an edge tile the canvas without gaps or double work.
#
@inline function _raster_triangle!(F::Matrix{Float64},
                                   x1::Float64, y1::Float64, v1::Float64,
                                   x2::Float64, y2::Float64, v2::Float64,
                                   x3::Float64, y3::Float64, v3::Float64,
                                   λ0::Float64, φ0::Float64,
                                   dλ::Float64, dφ::Float64,
                                   nlon::Int, nlat::Int)

    den = (y2 - y3)*(x1 - x3) + (x3 - x2)*(y1 - y3)
    #
    # A triangle of zero area in the (λ, φ) plane. Two of its nodes coincide, or
    # all three are collinear — which is what a triangle touching a POLE becomes,
    # since the pole is a whole line of the map. It paints nothing, and the
    # pixels it would have covered are the ones _fill_raster_gaps! picks up.
    #
    abs(den) < 1.0e-30 && return nothing
    iden = 1.0/den

    xmin = min(x1, x2, x3); xmax = max(x1, x2, x3)
    ymin = min(y1, y2, y3); ymax = max(y1, y2, y3)

    imin = max(1,    ceil(Int,  (xmin - λ0)/dλ + 0.5))
    imax = min(nlon, floor(Int, (xmax - λ0)/dλ + 0.5))
    jmin = max(1,    ceil(Int,  (ymin - φ0)/dφ + 0.5))
    jmax = min(nlat, floor(Int, (ymax - φ0)/dφ + 0.5))

    #
    # -1e-9 rather than 0: a pixel centre sitting exactly on a shared edge is
    # inside BOTH triangles in exact arithmetic and can be inside NEITHER in
    # floating point, which would leave a one-pixel crack along every element
    # edge. The overlap this tolerance creates is harmless — both triangles
    # agree on the edge.
    #
    tol = -1.0e-9

    @inbounds for i = imin:imax
        x = λ0 + (i - 0.5)*dλ
        for j = jmin:jmax
            y = φ0 + (j - 0.5)*dφ

            b1 = ((y2 - y3)*(x - x3) + (x3 - x2)*(y - y3))*iden
            b1 < tol && continue
            b2 = ((y3 - y1)*(x - x3) + (x1 - x3)*(y - y3))*iden
            b2 < tol && continue
            b3 = 1.0 - b1 - b2
            b3 < tol && continue

            F[i,j] = b1*v1 + b2*v2 + b3*v3
        end
    end
    return nothing
end


#
# Pixels no triangle covered — the polar caps (see the header) — take the value
# of the nearest node, measured as the CHORD in 3-D. The sphere has no polar
# singularity; only its map does, so the search is done off the map.
#
function _fill_raster_gaps!(F::Matrix{Float64}, λ::Vector{Float64}, φ::Vector{Float64},
                            f::AbstractVector, lon::AbstractVector, lat::AbstractVector,
                            npoin::Int, nfill_sample::Int)

    gaps = count(isnan, F)
    gaps == 0 && return F

    # Subsample: the search is brute force, and the caps are a few hundred
    # pixels of a 260 000-pixel canvas.
    stride = max(1, cld(npoin, max(nfill_sample, 1)))
    idx    = collect(1:stride:npoin)
    ns     = length(idx)

    sx = Vector{Float64}(undef, ns)
    sy = Vector{Float64}(undef, ns)
    sz = Vector{Float64}(undef, ns)
    sv = Vector{Float64}(undef, ns)
    @inbounds for k = 1:ns
        ip     = idx[k]
        sφ, cφ = sincos(Float64(lat[ip]))
        sλ, cλ = sincos(Float64(lon[ip]))
        sx[k]  = cφ*cλ
        sy[k]  = cφ*sλ
        sz[k]  = sφ
        sv[k]  = Float64(f[ip])
    end

    deg2rad_ = π/180.0
    @inbounds for j = 1:size(F,2)
        sφp, cφp = sincos(φ[j]*deg2rad_)
        for i = 1:size(F,1)
            isnan(F[i,j]) || continue
            sλp, cλp = sincos(λ[i]*deg2rad_)
            px, py, pz = cφp*cλp, cφp*sλp, sφp
            dbest = Inf; kbest = 1
            for k = 1:ns
                d = (px - sx[k])^2 + (py - sy[k])^2 + (pz - sz[k])^2
                if d < dbest
                    dbest = d; kbest = k
                end
            end
            F[i,j] = sv[kbest]
        end
    end
    return F
end
