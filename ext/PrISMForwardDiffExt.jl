module PrISMForwardDiffExt

using PrISM, ForwardDiff
import PrISM: find_c, get_c!, _dltar_c

function get_val_model(m::RWModel)
    RWModel(ForwardDiff.value.(m.m), ForwardDiff.value.(m.h),
        ForwardDiff.value.(m.ρ), ForwardDiff.value.(m.vp))
end

function get_val_model(m::LWModel)
    LWModel(ForwardDiff.value.(m.m), ForwardDiff.value.(m.h), ForwardDiff.value.(m.ρ))
end

function get_c!(resp_::AbstractArray{<:ForwardDiff.Dual{T, V, N}},
        t, m, mode, dc, c_start, c_high) where {T, V, N}
    m_val = get_val_model(m)
    for i in eachindex(t) # this can be parallelized
        ω = 2π / t[i]

        # c_low = copy(c_start)
        c_high_each = c_start

        for im in 1:(mode + 1)
            c = find_c(_dltar_c, ForwardDiff.value(c_high_each),
                ForwardDiff.value(c_high), ω, dc, m_val)
            c_high_each = ForwardDiff.value(c) + dc
        end
        c = c_high_each - dc
        fₓ = ForwardDiff.value(ForwardDiff.derivative(c_ -> _dltar_c(c_, ω, m_val), c))
        fₚ = ForwardDiff.partials(_dltar_c(c, ω, m))
        resp_[i] = ForwardDiff.Dual{T, V, N}(c, -fₚ / fₓ)
    end
    nothing
end

end
