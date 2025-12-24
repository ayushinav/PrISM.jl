"""
    mutable struct DCModelDistribution{T1<: Union{Distribution, AbstractArray}, T2<: Union{Distribution, AbstractArray}}
        m::T1
        h::T2
    end

creates a placeholder to store the `Distributions.jl` samplers for a priori
"""
mutable struct DCModelDistribution{
    T1 <: Union{Distribution, AbstractArray}, T2 <: Union{Distribution, AbstractArray}} <:
               AbstractGeophyModelDistribution{T1, T2}
    m::T1
    h::T2
end

"""
    struct DCResponseDistribution{T1<: Union{Function, Nothing}, T2<: Union{Function, Nothing}}
        ρₐ::T1
        ϕ::T2
    end

create a placeholder to store functions to obtain `Distributions.jl` samplers for the likelihood function
"""
struct DCResponseDistribution{T1 <: Union{Function, Nothing}} <:
       AbstractGeophyResponseDistribution
    ρₐ::T1
end
