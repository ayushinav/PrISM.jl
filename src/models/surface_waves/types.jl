mutable struct SurfaceWaveResponse{T} <: AbstractGeophyResponse
    c::T
end

mutable struct LWModel{
    T1 <: AbstractArray{<:Any}, T2 <: AbstractArray{<:Any}, T3 <: AbstractArray{<:Any}} <:
               AbstractGeophyModel{T1, T2}
    m::T1
    h::T2
    ρ::T3
end

mutable struct RWModel{T1 <: AbstractArray{<:Any}, T2 <: AbstractArray{<:Any},
    T3 <: AbstractArray{<:Any}, T4 <: AbstractArray{<:Any}} <: AbstractGeophyModel{T1, T2}
    m::T1
    h::T2
    ρ::T3
    vp::T4
end

RWModel(m, h, ρ, γ::T) where {T <: Number} = RWModel(m, h, ρ, m .* γ)

function SubsurfaceCore.default_params(::Type{T}) where {T <: Union{LWModel, RWModel}}
    (; mode=0, dc=0.01, dt=0.001, type=:phase)
end
const default_params_surface_waves = default_params(RWModel) # same for LWModel
