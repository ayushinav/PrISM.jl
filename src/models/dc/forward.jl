"""
get T for one wavelength
"""
function get_T(m, h, k)
    
    T = m[end]
    
    j = length(h)-1
    while j>=1
        T = (T + m[j]*tanh(k * h[j]))/ (1 + T * m[j]*tanh(k * h[j]))
    end

    return T

end


"""
we can have relative r_locs 
src_a_loc : 1 element
src_b_loc : 1 element
rec_m_loc : n, element
rec_n_loc : n, element

srcs = [src_a_loc, src_b_loc]
recs = hcat(rec_m_loc ; rec_n_loc) 2 x n
"""
function fwd(m, h, srcs, recs, hankel_filter)
    r_a = abs.(recs[:] .- srcs[1])
    r_b = abs.(recs[:] .- srcs[2])
    r_min, r_max = extrema(r_a..., r_b...)

    log_λ_min = log(first(hankel_filter.base)/r_max)
    log_λ_max = log(hankel_filter.base[end]/r_min)

    n_hankel = length(hankel_filter.base)

    d_log_λ = log(hankel_filter.base[1]) - log(hankel_filter.base[0])
    n_points = Int((log_λ_max ÷ d_log_λ) + 1) # TODO : for λ_max close λ_min
    log_λ_end = n_points * d_log_λ

    λs = exp.(range(log_λ_min, log_λ_end, length = n_points))
    nr = n_points + n_hankel
    log_r_end = nr * d_log_λ

    rs = r_min * exp.(range(0, log_r_end, length = nr))

    T = [get_T(m, h, λ) for λ in λs] # TODO

    appres = zeros(eltype(m), size(recs, 2))
    V_r = zero(appres)

    # V(r) = ∫ T(λ) J₀(λ r) dλ ≈ ∑ T_i * J₀_i
    # ρₐ = 
    for ir in eachindex(V_r)
        Vr[ir] = hankel_filter.J₀ ⋅ T[ir:ir+n_hankel-1]
    end

    
    # for ir in axes(recs, 2)
    #     Va = l;
    # end



end
