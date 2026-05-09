module PrISMForwardDiffExt

using PrISM, ForwardDiff
import PrISM: find_c_

function _ift_find_c(f, c1_val, c2_val, ω, m_val, m_dual,
        ::Type{ForwardDiff.Dual{T, V, N}}) where {T, V, N}
    c = find_c_(f, c1_val, c2_val, ω, m_val)
    fₓ = ForwardDiff.value(ForwardDiff.derivative(c_ -> f(c_, ω, m_val), c))
    fₚ = ForwardDiff.partials(f(c, ω, m_dual))
    return ForwardDiff.Dual{T, V, N}(c, -fₚ / fₓ)
end

# RWModel 

# all fields (rarely useful)
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
    _ift_find_c(f, ForwardDiff.value(c1), ForwardDiff.value(c2),
        ω, m_val, m, ForwardDiff.Dual{T, V, N})
end

# ∂m, ∂h
function find_c_(f,
        c1,
        c2,
        ω,
        m::RWModel{<:AbstractArray{<:ForwardDiff.Dual{T, V, N}},
            <:AbstractArray{<:ForwardDiff.Dual{T, V, N}},
            <:AbstractArray, <:AbstractArray}) where {T, V, N}
    m_val = RWModel(ForwardDiff.value.(m.m), ForwardDiff.value.(m.h), m.ρ, m.vp)
    _ift_find_c(f, ForwardDiff.value(c1), ForwardDiff.value(c2),
        ω, m_val, m, ForwardDiff.Dual{T, V, N})
end

# ∂m
function find_c_(f,
        c1,
        c2,
        ω,
        m::RWModel{<:AbstractArray{<:ForwardDiff.Dual{T, V, N}},
            <:AbstractArray, <:AbstractArray, <:AbstractArray}) where {T, V, N}
    m_val = RWModel(ForwardDiff.value.(m.m), m.h, m.ρ, m.vp)
    _ift_find_c(f, ForwardDiff.value(c1), ForwardDiff.value(c2),
        ω, m_val, m, ForwardDiff.Dual{T, V, N})
end

# ∂h
function find_c_(f,
        c1,
        c2,
        ω,
        m::RWModel{<:AbstractArray, <:AbstractArray{<:ForwardDiff.Dual{T, V, N}},
            <:AbstractArray, <:AbstractArray}) where {T, V, N}
    m_val = RWModel(m.m, ForwardDiff.value.(m.h), m.ρ, m.vp)
    _ift_find_c(f, ForwardDiff.value(c1), ForwardDiff.value(c2),
        ω, m_val, m, ForwardDiff.Dual{T, V, N})
end

# ∂ρ
function find_c_(f,
        c1,
        c2,
        ω,
        m::RWModel{<:AbstractArray, <:AbstractArray,
            <:AbstractArray{<:ForwardDiff.Dual{T, V, N}}, <:AbstractArray}) where {T, V, N}
    m_val = RWModel(m.m, m.h, ForwardDiff.value.(m.ρ), m.vp)
    _ift_find_c(f, ForwardDiff.value(c1), ForwardDiff.value(c2),
        ω, m_val, m, ForwardDiff.Dual{T, V, N})
end

# ∂vp
function find_c_(f,
        c1,
        c2,
        ω,
        m::RWModel{<:AbstractArray, <:AbstractArray, <:AbstractArray,
            <:AbstractArray{<:ForwardDiff.Dual{T, V, N}}}) where {T, V, N}
    m_val = RWModel(m.m, m.h, m.ρ, ForwardDiff.value.(m.vp))
    _ift_find_c(f, ForwardDiff.value(c1), ForwardDiff.value(c2),
        ω, m_val, m, ForwardDiff.Dual{T, V, N})
end

# LWModel 

# all fields (rarely useful)
function find_c_(f,
        c1,
        c2,
        ω,
        m::LWModel{<:AbstractArray{<:ForwardDiff.Dual{T, V, N}},
            <:AbstractArray{<:ForwardDiff.Dual{T, V, N}},
            <:AbstractArray{<:ForwardDiff.Dual{T, V, N}}}) where {T, V, N}
    m_val = LWModel(
        ForwardDiff.value.(m.m), ForwardDiff.value.(m.h), ForwardDiff.value.(m.ρ))
    _ift_find_c(f, ForwardDiff.value(c1), ForwardDiff.value(c2),
        ω, m_val, m, ForwardDiff.Dual{T, V, N})
end

# ∂m, ∂h
function find_c_(f,
        c1,
        c2,
        ω,
        m::LWModel{<:AbstractArray{<:ForwardDiff.Dual{T, V, N}},
            <:AbstractArray{<:ForwardDiff.Dual{T, V, N}}, <:AbstractArray}) where {T, V, N}
    m_val = LWModel(ForwardDiff.value.(m.m), ForwardDiff.value.(m.h), m.ρ)
    _ift_find_c(f, ForwardDiff.value(c1), ForwardDiff.value(c2),
        ω, m_val, m, ForwardDiff.Dual{T, V, N})
end

# ∂m
function find_c_(f,
        c1,
        c2,
        ω,
        m::LWModel{
            <:AbstractArray{<:ForwardDiff.Dual{T, V, N}}, <:AbstractArray, <:AbstractArray}) where {
        T, V, N}
    m_val = LWModel(ForwardDiff.value.(m.m), m.h, m.ρ)
    _ift_find_c(f, ForwardDiff.value(c1), ForwardDiff.value(c2),
        ω, m_val, m, ForwardDiff.Dual{T, V, N})
end

# ∂h
function find_c_(f,
        c1,
        c2,
        ω,
        m::LWModel{
            <:AbstractArray, <:AbstractArray{<:ForwardDiff.Dual{T, V, N}}, <:AbstractArray}) where {
        T, V, N}
    m_val = LWModel(m.m, ForwardDiff.value.(m.h), m.ρ)
    _ift_find_c(f, ForwardDiff.value(c1), ForwardDiff.value(c2),
        ω, m_val, m, ForwardDiff.Dual{T, V, N})
end

# ∂ρ
function find_c_(f,
        c1,
        c2,
        ω,
        m::LWModel{
            <:AbstractArray, <:AbstractArray, <:AbstractArray{<:ForwardDiff.Dual{T, V, N}}}) where {
        T, V, N}
    m_val = LWModel(m.m, m.h, ForwardDiff.value.(m.ρ))
    _ift_find_c(f, ForwardDiff.value(c1), ForwardDiff.value(c2),
        ω, m_val, m, ForwardDiff.Dual{T, V, N})
end

end
