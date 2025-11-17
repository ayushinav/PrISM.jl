"""
    mutable struct LWModelDistribution{T1<: Union{Distribution, AbstractArray}, T2<: Union{Distribution, AbstractArray}} # where T1,T2 
        m::T1
        h::T2
        ρ::T3
    end

creates a placeholder to store the `Distributions.jl` samplers for a priori
"""
mutable struct LWModelDistribution{
    T1 <: Union{Distribution, AbstractArray}, T2 <: Union{Distribution, AbstractArray},
    T3 <: Union{Distribution, AbstractArray}} <: AbstractGeophyModelDistribution{T1, T2}
    m::T1
    h::T2
    ρ::T3
end

"""
    mutable struct RWModelDistribution{T1<: Union{Distribution, AbstractArray}, T2<: Union{Distribution, AbstractArray}} # where T1,T2 
        m::T1
        h::T2
        ρ::T3
        vp:T4
    end

creates a placeholder to store the `Distributions.jl` samplers for a priori
"""
mutable struct RWModelDistribution{
    T1 <: Union{Distribution, AbstractArray}, T2 <: Union{Distribution, AbstractArray},
    T3 <: Union{Distribution, AbstractArray}, T4 <: Union{Distribution, AbstractArray}} <:
               AbstractGeophyModelDistribution{T1, T2}
    m::T1
    h::T2
    ρ::T3
    vp::T4
end

"""
    struct SurfaceWaveResponseDistribution{T1<: Union{Function, Nothing}}
        c::T1
    end

creates a placeholder to store the `Distributions.jl` samplers for a priori
"""
struct SurfaceWaveResponseDistribution{T1 <: Union{Function, Nothing}} <:
       AbstractGeophyResponseDistribution
    c::T1
end
