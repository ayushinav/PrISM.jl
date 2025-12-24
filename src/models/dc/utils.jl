"""
    hankel_transform_and_interpolation(r_max, r_min, fn, hankel_filter)

performs Hankel transform to evaluate the integral of the form:
```math
∫ fn(λ) J₀ (λr) dλ
```
and returns an interpolating function for any r ∈ [r_min, r_max]

### Arguments

 - `r_max` : maximum r 
 - `r_min` : minimum r
 - `fn` : kernel function
 - `hankel_filter` : hankel filter to be used, contains base, J₀ and J₁

returns an interpolated function that outputs the integral using Hankel transform at any r
"""
function hankel_transform_and_interpolation(r_max, r_min, fn, hankel_filter)
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
