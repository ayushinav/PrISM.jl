mutable struct multi_phase_modelType2{T1, T2, T3, T4, T5, T6, T7, T8, M}
    m1::Type{T1}
    m2::Type{T2}
    m3::Type{T3}
    m4::Type{T4}
    m5::Type{T5}
    m6::Type{T6}
    m7::Type{T7}
    m8::Type{T8}
    mix::M
end

multi_phase_modelType2(m1) = m1
multi_phase_modelType2(m1, m::MT.phase_mixing) = m1

check_fn(m1) = m1

for i in 2:7
    args = [Symbol("m$k") for k in 1:i]
    last_args = :(m::MT.phase_mixing)
    expr_lhs = Expr(:call, :multi_phase_modelType2, args..., last_args)

    args2 = [Nothing for k in i+1:8]
    expr_rhs = Expr(:call, :multi_phase_modelType2, args..., args2..., last_args)

    # expr = Expr(:(=), expr_lhs, expr_rhs)
    expr = Expr(:function, expr_lhs, expr_rhs)
    @show expr
    eval(expr)
end


@inferred multi_phase_modelType2(SEO3, Ni2011, Gaillard2008, HS1962_plus())

mutable struct multi_phase_model3{T, T1, T2, T3, T4, T5, T6, T7, T8, M}
    ϕ::T
    m1::T1
    m2::T2
    m3::T3
    m4::T4
    m5::T5
    m6::T6
    m7::T7
    m8::T8
    mix::M
end


function (model::multi_phase_modelType2)(ps::NamedTuple)
    pnames = propertynames(model)
    mix = getfield(model, pnames[end])
    # ϕ = getfield(ps, :ϕ)

    # make ϕ according to different number of phases
    # mvec = [getfield(model, Symbol("m$i")) for i in 1:8]
    ϕ_vec = rearrange_ϕ2(ps.ϕ, model)

    # pnames = pnames[2:3]
    v1 = from_nt(getproperty(model, pnames[1]), ps)
    v2 = from_nt(getproperty(model, pnames[2]), ps)
    v3 = from_nt(getproperty(model, pnames[3]), ps)
    v4 = from_nt(getproperty(model, pnames[4]), ps)
    v5 = from_nt(getproperty(model, pnames[5]), ps)
    v6 = from_nt(getproperty(model, pnames[6]), ps)
    v7 = from_nt(getproperty(model, pnames[7]), ps)
    v8 = from_nt(getproperty(model, pnames[8]), ps)
    return multi_phase_model3(ϕ_vec, v1, v2, v3, v4, v5, v6, v7, v8, mix)
end

import MT:from_nt
function rearrange_ϕ2(ϕ, model::multi_phase_modelType2)
    @assert sum(ϕ)≤1 "Σϕᵢ = $(sum(ϕ)) should be ≤ 1."

    fnames = propertynames(model)[1:(end - 1)]
    fnames = filter(f -> getfield(model, f) !== Nothing, fnames)
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

# ====

n1 = multi_phase_modelType(SEO3, Ni2011, Zhang2012, MT.HS_minus_multi_phase())

T = collect(1000:100:1400) .+ 273.0f0
Ch2o_ol = 1.0f0
Ch2o_m = 1000.0f0
Ch2o_opx = 100.0f0
Cco2_m = 10.0f0
ϕ = [0.1f0, 0.2f0]

ps_nt = (; ϕ=ϕ, T=T, Ch2o_ol=Ch2o_ol, Ch2o_m=Ch2o_m, Cco2_m=Cco2_m, Ch2o_opx)

m1 = n1(ps_nt)


forward(m1, [])

# plus  : Float32[-3.2222052, -2.105757, -1.5998105, -1.103057, -0.5914401]
# minus :  Float32[2.5088632, 2.6440232, 2.7609413, 2.8631563, 2.9532573]
