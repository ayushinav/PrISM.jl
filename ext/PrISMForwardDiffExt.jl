module PrISMForwardDiffExt

using PrISM, ForwardDiff
import PrISM: find_c_

# const _DA = AbstractArray{<:ForwardDiff.Dual}

# @ForwardDiff_frule find_c_(f, c1, c2, ω, m::RWModel{<:AbstractArray{<:ForwardDiff.Dual}, <:AbstractArray, <:AbstractArray, <:AbstractArray})

# @ForwardDiff_frule find_c_(f, c1, c2, ω, m::RWModel{ <:AbstractArray, <:AbstractArray{<:ForwardDiff.Dual}, <:AbstractArray, <:AbstractArray})

# @ForwardDiff_frule find_c_(f, c1, c2, ω, m::RWModel{ <:AbstractArray, <:AbstractArray, <:AbstractArray{<:ForwardDiff.Dual}, <:AbstractArray})

# @ForwardDiff_frule find_c_(f, c1, c2, ω, m::RWModel{ <:AbstractArray, <:AbstractArray, <:AbstractArray, <:AbstractArray{<:ForwardDiff.Dual}})

function find_c_(f,
        c1,
        c2,
        ω,
        m::RWModel{<:AbstractArray{<:ForwardDiff.Dual{T, V, N}},
            <:AbstractArray{<:ForwardDiff.Dual{T, V, N}},
            <:AbstractArray{<:ForwardDiff.Dual{T, V, N}},
            <:AbstractArray{<:ForwardDiff.Dual{T, V, N}}}) where {T, V, N}
    m_val = RWModel(ForwardDiff.value.(m.m), ForwardDiff.value.(m.h),
        ForwardDiff.value.(m.ρ), ForwardDiff.value.(m.vp))
    c = find_c_(f, ForwardDiff.value(c1), ForwardDiff.value(c2), ω, m_val)

    # @show c
    # @show typeof(m_val)
    # @show m_val

    fₓ = ForwardDiff.value(ForwardDiff.derivative(c_ -> f(c_, ω, m_val), c))  # ∂f/∂c
    fₚ = ForwardDiff.partials(f(c, ω, m))                   # ∂f/∂m · Δm, free from one eval

    # @show "CUSTOM AD"
    # @show fₚ
    # @show fₓ
    # (-fₚ / fₓ).values

    # return ForwardDiff.Dual{T,V,N}(c, (-fₚ / fₓ).values...)

    return ForwardDiff.Dual{T, V, N}(c, -fₚ / fₓ)
end

function find_c_(f,
        c1,
        c2,
        ω,
        m::LWModel{<:AbstractArray{<:ForwardDiff.Dual{T, V, N}},
            <:AbstractArray{<:ForwardDiff.Dual{T, V, N}},
            <:AbstractArray{<:ForwardDiff.Dual{T, V, N}}}) where {T, V, N}
    m_val = LWModel(ForwardDiff.value.(m.m), ForwardDiff.value.(m.h), ForwardDiff.value.(m.ρ))
    c = find_c_(f, ForwardDiff.value(c1), ForwardDiff.value(c2), ω, m_val)

    # @show c
    # @show typeof(m_val)
    # @show m_val

    fₓ = ForwardDiff.value(ForwardDiff.derivative(c_ -> f(c_, ω, m_val), c))  # ∂f/∂c
    fₚ = ForwardDiff.partials(f(c, ω, m))                   # ∂f/∂m · Δm, free from one eval

    # @show "CUSTOM AD"
    # @show fₚ
    # @show fₓ
    # (-fₚ / fₓ).values

    # return ForwardDiff.Dual{T,V,N}(c, (-fₚ / fₓ).values...)

    return ForwardDiff.Dual{T, V, N}(c, -fₚ / fₓ)
end

end
