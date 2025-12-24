"""
get T for one wavelength
"""
function get_T(m, h, k)
    T = m[end]

    j = length(h)
    while j >= 1
        T = (T + m[j] * tanh(k * h[j])) / (1 + T * tanh(k * h[j]) / m[j])
        j -= 1
    end

    return T
end

function SubsurfaceCore.forward(m::Tm, locs::T3, params) where {Tm <: DCModel, T3}
    r_a = abs.(recs[:] .- srcs[1])
    r_b = abs.(recs[:] .- srcs[2])
    r_min, r_max = extrema([r_a..., r_b...])

    fn(k) = get_T(m, h, k)

    V_fn = hankel_transform_and_interpolation(r_max, r_min, fn, params.hankel_filter)

    r_am = abs.(recs[:, 1] .- srcs[1])
    V_am = V_fn.(log.(r_am))
    r_an = abs.(recs[:, 2] .- srcs[1])
    V_an = V_fn.(log.(r_an))

    r_bm = abs.(recs[:, 1] .- srcs[2])
    V_bm = V_fn.(log.(r_bm))
    r_bn = abs.(recs[:, 2] .- srcs[2])
    V_bn = V_fn.(log.(r_bn))

    ρₐ = (V_am - V_bm - V_an + V_bn) ./ (inv.(r_am) - inv.(r_bm) - inv.(r_an) + inv.(r_bn))
end
