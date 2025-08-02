using Pkg
Pkg.activate("pkgs/MT.jl/wip")

using MT

const global boltz_k = 8.617e-5;
const global charge_e = 1.602e-19;

# models

import MT: AbstractModel, AbstractResponse
import MT: AbstractModelDistribution, AbstractResponseDistribution
abstract type AbstractRockphyModel <: AbstractModel end
abstract type AbstractRockphyResponse <: AbstractResponse end
abstract type AbstractRockphyModelDistribution <: AbstractModelDistribution end
abstract type AbstractRockphyResponseDistribution <: AbstractResponseDistribution end


mutable struct Rockphycond{T} <: AbstractRockphyResponse
    σ::T
end

# TODO : use of @concrete types

mutable struct SEO3{F} <: AbstractRockphyModel
    T::F
end

mutable struct Ni2011{F1,F2} <: AbstractRockphyModel
    T::F1
    Ch2o_m::F2
end

# ====== mixing laws =====

mutable struct HS_plus{T1<:AbstractRockphyModel,T2} <: AbstractRockphyModel
    models::Vector{T1}
    ϕ::Vector{T2}

    # test that ∑ϕ = 1
end


# TODO : change these to f-5 and f-19

"""
T is in K
"""
function forward(m::SEO3)
    fO₂ = 10.0^(-24441.9 / m.T + 13.296)

    # sT 
    bfe = 5.06e24 * exp(-0.357 * inv(boltz_k * m.T))
    bmg = 4.58e26 * exp(-0.752 * inv(boltz_k * m.T))
    ufe = 12.2e-6 * exp(-1.05 * inv(boltz_k * m.T))
    umg = 2.72e-6 * exp(-1.09 * inv(boltz_k * m.T))
    concFe = bfe + 3.33e24 * exp(-0.02 * inv(boltz_k * m.T)) * fO₂^(1 / 6)
    concMg = bmg + 6.21e30 * exp(-1.83 * inv(boltz_k * m.T)) * fO₂^(1 / 6)
    σ = concFe * ufe * charge_e + 2.0 * concMg * umg * charge_e

    return Rockphycond(σ)

end


"""
T is in K
"""
function forward(m::Ni2011)

    Tcorr = 1146.8 # K
    D = 0.006 # unitless, Partition coefficient {ol/melt}

    ls = 2.172 - (860.82 - 204.46 * sqrt(m.Ch2o_m / 1e4)) * inv(m.T - Tcorr)
    σ = 10.0^ls

    return Rockphycond(σ)

end


# We need to define forward(m::HS_plus)

function forward(m::HS_plus)
    σs = forward.(m.models)

    σ_max = maximum(σs)

    σ_plus = inv(sum(m.ϕ .* inv.(σ_max .+ σs))) - σ_max

    return σ_plus
end


m1 = SEO3(1000)
m2 = Ni2011(1000, 0.01)

forward(m1)

T = 1000;
hs_model = HS_plus([SEO3(T), Ni2011(T, 0.1)], [0.4, 0.6])


forward(hs_model)
typeof(m2)

forward(hs_model.models[1])

c = (; zip([:a, :b], [100, 10])...)

m_list = [SEO3, Ni2011];
p = (; zip([:T, :Ch2o_m], [Normal(1000.0, 100.0), Uniform(0.0, 1.0)])...)

propertynames(c)


function get_property_names(m::HS_plus)
    var_list = []

    for i in eachindex(m.models)
        i_list = [propertynames(m.models[i])...]
        model_name = string(typeof(m.models[i]).name.name)

        push!(var_list, Symbol.(model_name .* "___" .* string.(i_list))...)

    end

    # [
    #     [propertynames(m.models[i])...] 
    # for i in eachindex(m.models)]
    return var_list

end

get_property_names(hs_model)




function make_HS_dists(m::Vector{T}, p::NamedTuple) where {T}

    # if we take a dictionary/named tuple as input, we can propogate the same variable name to all the later models unless mentioned explicitly mentioned otherwise
    # feature of two models having one T, and others having another T is not present at the moment?

    # For now, let's have all the models share the parameters with the same names

    var_list = [keys(p)...]









end






#=
if we had just Ni2011 and SEO3, we could have made a bigger struct with 2 parameters, and called the respective forward functions in the forward call. 
    Now, we have to make a forward function for each type. 

    Let's make a HS_plus struct as it is now. This gives us the list of parameters


    The easiest route would be to create a struct on the fly, and its forward function as well.
    Why do we need a struct
        1) forward call dispatch
        2) get parameters
        3) sample parameters and make the new model

    Now, we can overload to get the parameter list. We can then constuct the same HS_model by just plugging in the values.
    Forward is then easy.



=#

# overload alg_cache.apriori, ie, AbstractModelDistribution
# overload 



mutable struct HS_2
    T::Any
    ch2om::Any
    ϕ::Any
end


function forward(m::HS_2)
    m1 = SEO3(m.T)
    m2 = Ni2011(m.T, m.ch2om)

    σs = [m1, m2]

    σ_max = maximum(σs)

    σ_plus = inv(sum(m.ϕ .* inv.(σ_max .+ σs))) - σ_max

    return Rockphycond(σ_plus)

end


function fn(m; args...) # T = 1000, Ch2o_m = 0.1

    for im in m
        var_list = [fieldnames(im)...]
        im(args[var_list...]...)
    end

    # @show args, typeof(args)
    # @show args[:T]

end


fn([SEO3, Ni2011]; T = 1000, Ch2o_m = 0.1)


c = (; zip([:T, :Ch2o_m], [1000, 0.1])...)

ll = [propertynames(m2)...]

c[ll]

m_list = [];
for m in [m1, m2]
    var_list = [propertynames(m)...]
    push!(m_list, typeof(m)(c[var_list]...))
end


function fn(args...)
    @show args
end


fn(100, 10000)


mutable struct mixing_models5{T<:Vector{<:Any}}
    p::Any # vector of parameters 
    ϕ::Any # phase ratios
end


mixing_models4{[SEO3, Ni2011]}(c, [0.8, 0.2])



mutable struct check{T}
    m::T
end


c = check([100.0, 100.0, 10000])

typeof(c)
