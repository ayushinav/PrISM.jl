"""
    mutable struct modelDistribution{T1<: Union{Distribution, AbstractArray}, T2<: Union{Distribution, AbstractArray}} # where T1,T2 
        m::T1
        h::T2
    end

create a placeholder to store the `Distributions.jl` sampler for a priori
"""
mutable struct LWModelDistribution{
    T1 <: Union{Distribution, AbstractArray}, T2 <: Union{Distribution, AbstractArray}, T3 <: Union{Distribution, AbstractArray}} <:
               AbstractGeophyModelDistribution{T1, T2}
    m::T1
    h::T2
    ρ::T3
end

mutable struct RWModelDistribution{
    T1 <: Union{Distribution, AbstractArray}, T2 <: Union{Distribution, AbstractArray}, T3 <: Union{Distribution, AbstractArray}, T4 <: Union{Distribution, AbstractArray}} <:
               AbstractGeophyModelDistribution{T1, T2}
    m::T1
    h::T2
    ρ::T3
    vp::T4
end

"""
    struct responseDistribution{T1<: Union{Function, Nothing}, T2<: Union{Function, Nothing}} # where T1,T2
        ρₐ::T1
        ϕ::T2
    end

create a placeholder to store functions to obtain `Distributions.jl` samplers for the likelihood function
"""
struct SurfaceWaveResponseDistribution{
    T1 <: Union{Function, Nothing}, T2 <: Union{Function, Nothing}} <:
       AbstractGeophyResponseDistribution
    c::T1
end
