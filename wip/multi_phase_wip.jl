using Pkg
Pkg.activate(".")

using Distributions
import Distributions:Distribution
using MT

using Test

T = collect(1000:100:1400) .+ 273.0f0
Ch2o_ol = 1.0f0
Ch2o_m = 1000.0f0
Cco2_m = 10.0f0
ϕ = 0.1f0

m = two_phase_modelType(Yoshino2009, Sifre2014, HS1962_plus())

ps_nt = (; ϕ=ϕ, T=T, Ch2o_ol=Ch2o_ol, Ch2o_m=Ch2o_m, Cco2_m=Cco2_m)
model = m(ps_nt);

@inferred m(ps_nt)
@inferred forward(model, [])

mutable struct construct_multi_phase_model2{T<: Tuple, M}
    m_vec::T
    mix::M
end

mutable struct multi_phase_model{T1, T, M}
    ϕ::T1
    m_vec::T
    mix::M
end

@inferred construct_multi_phase_model2((SEO3, Ni2011), HS1962_minus())

m1 = construct_multi_phase_model2((SEO3, Ni2011), HS1962_minus())

function rearrange_ϕ(ϕ, model::construct_multi_phase_model2)
    @assert sum(ϕ)≤1 "Σϕᵢ = $(sum(ϕ)) should be ≤ 1."

    fnames = propertynames(model)[1:(end - 1)]
    fnames = model.m_vec #filter(f -> getfield(model, f) !== Nothing, fnames)
    msg = """
    Vol. frac of last component is defined automatically once the others are defined. 
    Make sure the length of porosity parameter `ϕ` is one less than the length of number of components.
    length of ϕ = $(length(ϕ))
    length of components = $(length(fnames))
    Check out the relevant documentation.
    """

    @assert length(ϕ)==(length(fnames) - 1) msg

    c = length(ϕ)
    ϕ_vec = zeros(eltype(ϕ), length(fnames))

    ϕ_vec[1:c] .= ϕ
    ϕ_vec[c + 1] = 1 - sum(ϕ)
    return ϕ_vec
end


rearrange_ϕ([0.1], m1)

function (model::construct_multi_phase_model2)(ps::NamedTuple)
    # pnames = propertynames(model)
    # mix = getfield(model, pnames[end])
    mix = model.mix
    # ϕ = getfield(ps, :ϕ)

    # make ϕ according to different number of phases
    ϕ_vec = rearrange_ϕ(ps.ϕ, model)

    # m1 = model.m_vec[1]

    # pnames = pnames[2:3]
    # m_vec = [MT.from_nt(m_i, ps) for m_i in model.m_vec]
    # m_vec = map(m1.m_vec) do x
    #     MT.from_nt(x, ps)
    # end

    m_vec = supertype(model.m_vec[1])[]
    for i in eachindex(model.m_vec)
        push!(m_vec, MT.from_nt(model.m_vec[i], ps), )
    end

    # m_vec = helper_make_nt2(model.m_vec, ps)
    @show length(m_vec)
    # m2 = MT.from_nt(m1, ps)
    # v1 = from_nt(getproperty(model, pnames[1]), ps)
    # v2 = from_nt(getproperty(model, pnames[2]), ps)
    # v3 = from_nt(getproperty(model, pnames[3]), ps)
    # v4 = from_nt(getproperty(model, pnames[4]), ps)
    # v5 = from_nt(getproperty(model, pnames[5]), ps)


    # return multi_phase_model(ϕ_vec, m_vec, mix)

    return map(x -> MT.from_nt(x, ps), model.m_vec.parameters)
end

function stable_creator(types_tuple::Type{T}, ps) where {T <: Tuple}
    # T.parameters gives us the tuple of concrete types, e.g., (Float32, Int64)
    # We can map the constructor for each type over this tuple of types.
    # The compiler can see everything and will unroll this map,
    # generating specialized code for each constructor call.
    return map(x -> MT.from_nt(x, ps), T.parameters)
end

# 1. Your input can still be a vector
my_types_vec = [Float32, Int64]
my_params = (1,)

# 2. Convert the vector of types into a tuple of types for dispatch
#    The `Tuple()` constructor and the splat operator `...` do this.
types_as_tuple = Tuple(my_types_vec)
typeof(types_as_tuple)

# 3. Pass the *type* of the tuple to our stable function
#    Note the use of `typeof()`
result = stable_creator(typeof(types_as_tuple), my_params)

[Float64, Int32]
types2 = Tuple([SEO3{Float64}, Ni2011{Float64, Float64}])


stable_creator(typeof(types2), ps_nt)
m1_ = m1(ps_nt);

@inferred m1(ps_nt)
@code_warntype m1(ps_nt)

map(m1.m_vec) do x
    MT.from_nt(x, ps)
end

m3 = two_phase_modelType(SEO3, Ni2011, HS1962_minus())

@generated function helper_make_nt2(m_tup::NTuple{N, K}, ps) where {N, K}
    
    # kv = :([:(MT.from_nt(im[$i], ps)) for i in 1:N])
    kv = [:(MT.from_nt(m_tup[$i], ps)) for i in 1:N]
    @show Expr(:block, :[$(kv...)])

    Expr(:block, :[$(kv...)])
end

@inferred 
helper_make_nt2(tp1, ps_nt);

tp1 = (SEO3, Ni2011)
typeof(tp1)

NTuple{2, UnionAll}


typeof((Type{SEO3}, Type{Ni2011}))






