using Pkg
Pkg.activate(".")

using MT
using Test
using Distributions
using Turing

n1 = multi_phase_modelType(SEO3, Sifre2014, Zhang2012, MT.HS_minus_multi_phase())

T = [1200f0] #collect(1000:100:1400) .+ 273.0f0
Ch2o_ol = 1.0f0
Ch2o_m = 1000.0f0
Cco2_m = 1000.0f0
Ch2o_opx = 100.0f0
Cco2_m = 10.0f0
ϕ = [0.1f0, 0.2f0]

ps_nt = (; ϕ=ϕ, T=T, Ch2o_ol=Ch2o_ol, Ch2o_m=Ch2o_m, Cco2_m=Cco2_m, Ch2o_opx)

m1 = n1(ps_nt)

pvec_ = [0.2, 0.3, 0.5]
cvec_ = Tuple([0.4, 1.0, 0.6])
@inferred MT.mix_models(cvec_, pvec_, HS_plus_multi_phase())

@inferred broadcast((sig) -> MT.mix_models(sig, pvec_, HS_plus_multi_phase()), cvec_...)

MT.forward(m1, [])


@inferred forward(m1, [])
@code_warntype forward(m1, [])


for i in 2:5
    fnames = [Symbol("m$j") for j in 1:i]
    expr = :(function mix_1(σ_tup::NTuple{$i,T}, phi, mix) where {T}
        broadcast((m1, m2, m3) -> MT.mix_models((m1, m2, m3), phi, mix), σ_tup...)
    end)
end

function mix_1(σ_tup::NTuple{3,T}, phi, mix) where {T}
    broadcast((m1, m2, m3) -> MT.mix_models((m1, m2, m3), phi, mix), σ_tup...)
end

expr = :(function mix_1(σ_tup::NTuple{3,T}, phi, mix) where {T}
        broadcast((m1, m2, m3) -> MT.mix_models((m1, m2, m3), phi, mix), σ_tup...)
end)

i_  = 3
fnames = [Symbol("m$j") for j in 1:i_]

eval(expr)

expr

expr.head

expr.args
expr.args[1]
expr.args[2]

expr.args[2].head

expr.args[2].args[1]
expr.args[2].args[2]
expr.args[2].args[3]

expr.args[2].args[3].head


expr.args[2].args[3].args
expr.args[2].args[3].args[1]
expr.args[2].args[3].args[2]
expr.args[2].args[3].args[3]



expr2 = Expr(
    :function,

)


@inferred mix_1(cvec_, pvec_, mix)

@generated function mix_model_help(σ_tup::NTuple{N,T}, phi::P, m::M) where {N,T,P,M}
    # Generate N unique symbols (e.g., :m1, :m2) for the anonymous function's arguments.
    fnames = [QuoteNode("m$i") for i in 1:N]
    σ_v = [:(σ_tup[$i]) for i in 1:N]

    # fnames = @macroexpand @nexprs
    # @show fnames
    # @show σ_v

    # Construct the expression that will be executed at runtime.
    expr = :(
        broadcast(
            # 1. Create a type-stable anonymous function with N arguments.
            # ($(fnames...)) -> MT.mix_models(($(fnames...),), phi, m),

            (($(fnames[1]), $(fnames[2]))) -> MT.mix_models($(fnames[1], fnames[2]), phi, m),

            # 2. At runtime, splat the input tuple `σ_tup` to feed its vectors to broadcast.
            $(σ_v[1], σ_v[2])

            # $[:(σ_tup[$i]) for i in 1:N]

            # [1., 2., 4.]
        )
    )

    # expr = :(sum(phi))

    # expr = :(function new_fn( $(fnames...) )
    #     # sum( QuoteNode($(fnames...)) )
    #     println("YES")
    # end)

    # expr = :(MT.mix_models(($(fnames...),), phi, m))
    # expr = :(MT.mix_models(σ_tup... , phi, m))
    # expr = :(MT.mix_models($(σ_v...) , phi, m))
    # Expr(:block, :[$(kv...)])

    @show expr

    return expr
end


mix_model_help(cvec_, pvec_, HS_plus_multi_phase())

fnames = [Symbol("m$i") for i in 1:3]

mix = HS_plus_multi_phase()

import MT:mix_models

expr_ = Expr(
    :call,
    [
        broadcast,
        [
            Expr(
                :->,
                [(fnames...)],
                [
                    Expr(:call,
                    [:mix_models,
                    ((fnames...)),
                    pvec_,
                    mix
                    ]
                )
                    # mix_models(((fnames...)), pvec_, mix)]
                ]
            )
        ],
        [
            cvec_
        ]
    ]


)


eval(expr_)

fnx_expr = quote
    # function mix_model_help2(σ_tup, phi, m)
        broadcast((sig...) -> MT.mix_models(sig, phi, m), σ_tup...)
    # end
end

fnx_expr.args[1]

fnx_expr.args[2].head

fnx_expr.args[2].args

fnx_expr.args[2].args[1]

fnx_expr.args[2].args[2]

fnx_expr.args[2].args[2].head
fnx_expr.args[2].args[2].args

fnx_expr.args[2].args[3]


# eval(fnx_expr)

@macroexpand Base.Cartesian.@ncall 5 MT.mix_models x
@macroexpand Base.Cartesian.@ncall(5, MT.mix_models, x)

import Base.Cartesian:@ncall
using Base.Cartesian


function help3(sig::Tuple{T1, T2, T3, T4, T5, T6, T7, T8}, phi::P, mix::M) where {
    T1, T2, T3, T4, T5, T6, T7, T8 <: Nothing, P, M}
    help3()
end

function check2(x::Tuple{T1, T2, T3}) where{T1, T2, T3}

    @show T1
end


check2(Tuple(cvec_))

expr_ = @macroexpand @nexprs 3 i -> cvec_[i]

k = 2

a1 = 10.
a2 = ones(1,10)
a3 = ones(1,1,10)

a_vec = [a1, a2, a3]

function check_fn(a::Tuple{T1, T2, T3}, b, c) where {T <: Nothing} 
    check_fn(a[1:end-1],b,c)
end

check_fn(a...) = +(a...)
# check_fn(a,b) = a + b

broadcast((x1, x2, x3) -> +(x1, x2, x3), a_vec...)


@macroexpand @ncall  MT.mix_models i->sig[1] 

Tuple([1,2,3,4])[1:3]


+(@nexprs k i -> pvec_[i])

func(x) = Tuple(x)
@macroexpand @ncall 3 func i -> cvec_[i]

@ncall 3 maximum i -> cvec_[i]

broadcast((sig...) -> mix_models(sig, model.ϕ, model.mix), σ_vec...)

func2 = x -> Tuple(x...)

@macroexpand @ncall 2 func i->cvec_[i]
@macroexpand @ncall 2 func2 i->cvec_[i]

# mix_model_help2(cvec_, pvec_, HS_plus_multi_phase())

mix_model_help(cvec_, pvec_, HS_plus_multi_phase())

@inferred forward(m1, [])
@code_warntype forward(m1, [])



:(broadcast(((m1, m2, m3)->begin
              #= /Users/asingh933/Desktop/pkgs/MT.jl/wip/mulit_phase_wip3.jl:41 =#
              MT.mix_models((m1, m2, m3), phi, m)
          end), σ_tup...))
fnames = [:m1, :m2, :m3]
expr = :(broadcast(((m1, m2, m3)-> MT.mix_models((m1, m2, m3), phi, m)), σ_tup...))

#


# plus  : Float32[-3.2222052, -2.105757, -1.5998105, -1.103057, -0.5914401]
# minus :  Float32[2.5088632, 2.6440232, 2.7609413, 2.8631563, 2.9532573]

m = multi_phase_modelType(SEO3, Sifre2014, Zhang2012, MT.HS_plus_multi_phase())

model = m(ps_nt)
resp = forward(model, [])
err_resp = RockphyCond(0.01 .* abs.(resp.σ))

m = MT.multi_phase_modelDistributionType(SEO3Distribution, Sifre2014Distribution, Zhang2012Distribution, MT.HS_plus_multi_phase())
ps_nt_dist = (;
    T=MvNormal([1000.0], [400.0;]),
    Ch2o_m=product_distribution([Uniform(50.0, 150.0)]),
    Cco2_m=[100.0], # ϕ=product_distribution([Uniform(0.01, 0.2)]),
    Ch2o_opx=product_distribution([Uniform(50.0, 150.0)]),
    ϕ = product_distribution([Uniform(0.0f0, 0.1f0), Uniform(0.0f0, 0.3f0)])
)

mdist = m(ps_nt_dist);
rdist = RockphyCondDistribution(MT.normal_dist)

m_cache = mcmc_cache(mdist, rdist, 10, Prior());
mcmc_chain_prior = stochastic_inverse(resp, err_resp, [], m_cache);

m_cache = mcmc_cache(mdist, rdist, 5, NUTS());
mcmc_chain_posterior = stochastic_inverse(resp, err_resp, [], m_cache);


using CairoMakie

f = Figure(size = (400, 1000))

labels = ["ϕ₁", "ϕ₂", "T (K)", "melt water conc. (ppm)", "opx water conc. (ppm)"]

for i in 1:5
    axi = Axis(f[i, 1]; title=labels[i])

    idx = mcmc_chain_prior.name_map[1][i]

    dat_prior = mcmc_chain_prior[idx][:]
    dat_posterior = mcmc_chain_posterior[idx][:]

    density!(axi, dat_prior; label="prior")
    density!(axi, dat_posterior; label="posterior")
end

f[1:5, 2] = Legend(f, f.content[end]; orientation=:vertical)
f


# model_list = get_model_list(mcmc_chain_posterior, mdist)

# ====================================================================


mix = MAL(0.5)
# mix = HS_plus_multi_phase()
# mix = GAL([0.4, 0.6])
# m_mix = multi_phase_modelType(SEO3, Sifre2014, Zhang2012, mix)
# m_mix = two_phase_modelType(SEO3, Sifre2014, mix)

m_mix = two_phase_modelType(SEO3, Sifre2014, mix)
m = multi_rp_modelType(typeof(m_mix), anharmonic, Nothing, Nothing)

# @inferred MT.multi_rp_modelType(typeof(m_mix), anharmonic, Nothing, Nothing)


T = [1200f0] #collect(1000:100:1400) .+ 273.0f0
Ch2o_ol = [1.0f0]
Ch2o_m = [1000.0f0]
Cco2_m = [1000.0f0]
Ch2o_opx = [100.0f0]
Cco2_m = [10.0f0]
# ϕ = [0.3f0] #[0.1f0, 0.2f0]
ϕ = [0.1f0, 0.2f0]
P = [4f0]
ρ = [3300f0]
m_MAL = 0.6f0
m_GAL = [0.4, 0.6]

ps_nt = (; ϕ=ϕ, T=T, Ch2o_ol=Ch2o_ol, Ch2o_m=Ch2o_m, Cco2_m=Cco2_m, Ch2o_opx, P, ρ, m_GAL, m_MAL)

m1 = m(ps_nt)

# todo : from_nt

# plus  : Float32[-3.2222052, -2.105757, -1.5998105, -1.103057, -0.5914401]
# minus :  Float32[2.5088632, 2.6440232, 2.7609413, 2.8631563, 2.9532573]

# m = multi_phase_modelType(SEO3, Sifre2014, Zhang2012, MT.HS_plus_multi_phase())
# m = two_phase_modelType(SEO3, Sifre2014, mix)

model = m(ps_nt)
resp = forward(model, [])
err_resp = multi_rp_response(RockphyCond(0.01 .* abs.(resp.cond.σ)), 
            RockphyElastic(
                0.01 .* abs.(resp.elastic.G),
                0.01 .* abs.(resp.elastic.G),
                0.01 .* abs.(resp.elastic.Vp),
                0.01 .* abs.(resp.elastic.Vs),
            ), Nothing, Nothing)

m_mixdist = MT.multi_phase_modelDistributionType(SEO3Distribution, Sifre2014Distribution, Zhang2012Distribution, mix)
# m_mixdist = MT.multi_phase_modelDistributionType(SEO3Distribution, Sifre2014Distribution, Zhang2012Distribution, GAL(product_distribution([Uniform(0.3f0, 0.7f0), Uniform(0.3f0, 0.7f0)])))

# m_mixdist = two_phase_modelDistributionType(SEO3Distribution, Sifre2014Distribution, MAL([0.3f0])) #[Uniform(0.3f0, 0.7f0)])))

# m_mixdist = two_phase_modelDistributionType(
#     SEO3Distribution, Sifre2014Distribution, HS1962_plus())
m = multi_rp_modelDistributionType(
    typeof(m_mixdist), anharmonicDistribution, Nothing, Nothing)
ps_nt_dist = (;
    T=product_distribution([Uniform(1000., 1400.)]),
    Ch2o_m=product_distribution([Uniform(50.0, 150.0)]),
    Cco2_m=[100.0], # ϕ=product_distribution([Uniform(0.01, 0.2)]),
    Ch2o_opx=product_distribution([Uniform(50.0, 150.0)]),
    # ϕ = product_distribution([Uniform(0.0f0, 0.1f0)]),
    ϕ = product_distribution([Uniform(0.0f0, 0.1f0), Uniform(0.0f0, 0.3f0)]),
    P = [4f0], ρ = [3300f0],
    m_MAL = [0.3f0], #product_distribution([Uniform(0.3f0, 0.7f0)])
    m_GAL = [0.4, 0.6]
    # m_GAL = product_distribution([Uniform(0.3f0, 0.7f0), Uniform(0.3f0, 0.7f0)])
)

mdist = m(ps_nt_dist);
@show typeof(mdist)
rdist = MT.multi_rp_responseDistribution(
    RockphyCondDistribution(normal_dist),
    RockphyElasticDistribution(normal_dist, normal_dist, normal_dist, normal_dist),
    Nothing, Nothing
)

m_cache = mcmc_cache(mdist, rdist, 50, Prior());
mcmc_chain_prior = stochastic_inverse(resp, err_resp, [], m_cache);

m_cache = mcmc_cache(mdist, rdist, 20, NUTS());
mcmc_chain_posterior = stochastic_inverse(resp, err_resp, [], m_cache);


MT.to_nt(mix)

k_ = MAL{Uniform{Float32}}


function check(m::Type{T}) where T
    @show T
end

check(m_mix)

m3 = MT.from_nt(Zhang2012, ps_nt)
@inferred forward(m3, [])
@code_warntype forward(m3, [])