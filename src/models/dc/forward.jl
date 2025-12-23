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

"""
r_max : maximum r
r_min : minimum r
fn : kernel function
hankel_filter : filter to be used, contains base and J₀

returns an interpolated function that outputs the integral using Hankel transform at any r
"""
function do_Hankel_stuff(r_max, r_min, fn, hankel_filter)
    base_ = hankel_filter.base
    j0 = hankel_filter.J₀
    n_filter = length(base_)
    f = log(base_[2]) - log(base_[1])

    n_r = Int(ceil((log(r_max) - log(r_min)) / f)) + 1
    log_r_end = log(r_min) + (n_r - 1) * f # a + (n-1)d on log-scale

    rs = exp.(range(log_r_end, log(r_min); length=n_r))
    nλ = n_filter + n_r - 1

    log_λ_max = log(base_[end] / r_min)
    log_λ_start = log_λ_max - (nλ - 1) * f
    λs = exp.(range(log_λ_start, log_λ_max; length=nλ))

    T_at_λs = fn.(λs)
    T_at_rs = zeros(n_r)

    for ir in 1:n_r
        T_at_rs[ir] = inv(rs[ir]) .* (j0 ⋅ T_at_λs[ir:(ir + n_filter - 1)])
    end

    # make interpolating function and return

    f_spline = cubic_spline_interpolation(
        range(log(r_min), log_r_end; length=n_r), reverse(T_at_rs))

    return f_spline
end
"""
we can have relative r_locs
src_a_loc : 1 element
src_b_loc : 1 element
rec_m_loc : n, element
rec_n_loc : n, element

srcs = [src_a_loc, src_b_loc]
recs = hcat(rec_m_loc ; rec_n_loc) n x 2
"""
function fwd(m, h, srcs, recs, hankel_filter)
    r_a = abs.(recs[:] .- srcs[1])
    r_b = abs.(recs[:] .- srcs[2])
    r_min, r_max = extrema([r_a..., r_b...])

    fn(k) = get_T(m, h, k)

    V_fn = do_Hankel_stuff(r_max, r_min, fn, hankel_filter)

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
