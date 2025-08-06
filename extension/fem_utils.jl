"""
xrange= vector to be updated, of the same length as Δx
Δx= interval to the left of the grid point. The first is 0, making it the same length as the number of elements.
"""
function build_grid!(xrange, Δx)
    for i = 2:length(Δx)
        xrange[i] = xrange[i-1] + Δx[i]
    end
end

"""
`get_grid_for_layer(L; n= 16)` returns a set of grid points in a layer such that the separation between the points increases arithmatically towards the middle and then decreases after that with the same slope where `L` is the thickness of the layer constructing 2`n`-1 grid points.
"""
function get_grid_for_layer(L; n = 16)
    d = L / (n - 1) / n
    arr = zeros(n)
    an = 0
    for i = 1:n
        arr[i] = an
        an += d
    end
    b = arr[2:end]
    return [b..., reverse(b)...]
end

"""
`make_dz(ρ, h, ω, mzp; n= 16)` discretizes the resistivity distribution `ρ` with thickness `h` where `mzp` is used to get the maximum extent of the grid space for the frequency `ω` constructing 2`n`-1 grid points.
"""
function make_dz(ρ, h, ω, mzp, dzp; n = 16)
    n2 = 2n - 2
    max_depth = 500 * sqrt(maximum(ρ[2:end]) * 2π / ω) * mzp
    ngrid = (n2) * (length(h) + 1)
    # arr= zeros(n2*length(h)+n2);
    arr = Float64[]
    for i in eachindex(h)
        min_dz = 500 * sqrt(minimum(ρ[i]) * 2π / ω) * dzp
        n = max(n, Int(ceil(h[i] / min_dz)))
        # arr[n2*(i-1)+1:n2*(i)].= get_grid_for_layer(h[i], n=n);
        arr_ = get_grid_for_layer(h[i], n = n)
        push!(arr, arr_...)
    end
    # arr[end- n2+1: end].= get_grid_for_layer(max_depth- sum(h), n=n);
    arr_ = get_grid_for_layer(max_depth - sum(h), n = n)
    push!(arr, arr_...)

    return arr
end

function make_dz(ρ, h, ω::AbstractVector{T}, mzp, dzp; n = 16) where {T}
    n2 = 2n - 2
    max_depth = 500 * sqrt(maximum(ρ[2:end]) * 2π / minimum(ω)) * mzp
    # min_dz = 500*sqrt(minimum(ρ[2:end])*2π/ω)*dzp;
    ngrid = (n2) * (length(h) + 1)
    # arr= zeros(n2*length(h)+n2);
    arr = Float64[]
    for i in eachindex(h)
        min_dz = 500 * sqrt(minimum(ρ[i]) * 2π / maximum(ω)) * dzp
        n = max(n, Int(ceil(h[i] / min_dz)))
        # arr[n2*(i-1)+1:n2*(i)].= get_grid_for_layer(h[i], n=n);
        arr_ = get_grid_for_layer(h[i], n = n)
        push!(arr, arr_...)
    end
    # arr[end- n2+1: end].= get_grid_for_layer(max_depth- sum(h), n=n);
    arr_ = get_grid_for_layer(max_depth - sum(h), n = n)
    push!(arr, arr_...)

    return arr
end

function get_grids(ρ_, h_, ω, mzp, dzp; n = 32, air_thickness = 100)
    # @show extrema(ρ_)
    ρ = [1e12, ρ_...]
    # ρ = ρ_ #[1e12, ρ_...]
    h = [air_thickness, h_...]
    # h = h_ #[air_thickness, h_...]

    dz = make_dz(ρ, h, ω, mzp, dzp, n = n) # includes padding layer
    # @show dz
    zgrid = zero(dz)

    build_grid!(zgrid, dz)
    z_points = (zgrid[2:end] .+ zgrid[1:end-1]) ./ 2
    # @show extrema(z_points), z_points[1:5]
    # @show extrema(zgrid), zgrid[1:5]

    ρgrid = get_ρ_at_z([ρ..., h...], z_points)
    zgrid .= zgrid .- air_thickness

    return zgrid, ρgrid
end
