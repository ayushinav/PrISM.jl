mutable struct s1{T,K}
    m_vec::T
    mix::K
end
mix_ = HS1962_plus()
mixt_ = HS1962_plus
so1 = s1((SEO3, Ni2011), mix_)

s1{(SEO3, Ni2011), mixt_}
s1{(Float64, Float32), mixt_}

mutable struct s2{T <: Tuple,K}
    m_vec::T
    mix::K
end

abstract type nphasetype{T <: Tuple, K} end
mutable struct nphasetype1{T <: Tuple, K} 
end

mutable struct nphasetype3{T <: NamedTuple, K} 
end

function nphasetype3(tup::T, K) where {T <: Tuple}
    nt = NamedTuple{
        Tuple(Symbol("m$i") for i in eachindex(tup)),
        T}

    # ntt = @NamedTuple{}
    nphasetype3{nt, K}
end

mutable struct nphasetype6{T <: NamedTuple, K} 
    mvec::T
    mix::K
end


function nphasetype6(tup::T, mix) where {T <: Tuple}
    nt = (;zip(
        Tuple(Symbol("m$i") for i in eachindex(tup)),
        tup)
        ...)

    # ntt = @NamedTuple{}
    nphasetype6(nt, mix)
end

mutable struct nphasemodel{T <: Tuple, K,V}
    m_vec::T
    mix::K
    ϕ::V
end

mutable struct nphasemodel1{T, K, V}
    m_vec::T
    mix::K
    ϕ::V
end

n1 = nphasetype1{Tuple{Type{SEO3}, Type{Ni2011}}, mixt_}()
n1 = nphasetype1{Tuple{SEO3, Ni2011}, mixt_}()
# n1 = nphasetype2((SEO3, Ni2011), mix_)

function (m::nphasetype1{T,K})(ps) where {T,K}
    
    m_vec = map(T.types) do x
        MT.from_nt(x, ps)
    end
    # ϕ_vec = rearrange_ϕ(ps.ϕ, model)

    # stype = supertype(T.types[end])
    # nt_type = NTuple{length(T.types), stype}

    # m_vec = ntuple(
    #     i -> MT.from_nt(T.types[i], ps_nt), length(T.types)
    # )

    # @show T.parameters
    # @show T.parameters[1]

    # npm = MT.from_nt(T.parameters[1].parameters[1], ps_nt)
    npm = nphasemodel1(m_vec, K, ϕ)

    return npm
end

@inferred n1(ps_nt)
mo1 = n1(ps_nt)

function (m::nphasetype6{T,K})(ps) where {T,K}
    
    # m_vec = map(T.types) do x
    #     MT.from_nt(x, ps)
    # end
    # ϕ_vec = rearrange_ϕ(ps.ϕ, model)

    # stype = supertype(T.types[end])
    # nt_type = NTuple{length(T.types), stype}

    # m_vec = ntuple(
    #     i -> MT.from_nt(T.types[i], ps_nt), length(T.types)
    # )

    # @show T.parameters
    # @show T.parameters[1]
    @show T
    ktype = T.parameters[2].parameters[1]
    @show ktype

    # npm = MT.from_nt(ktype, ps_nt)
    npm = ktype

    #nphasemodel(Tuple(m_vec), K, ϕ)

    return npm
end

# MT.from_nt(Type{SEO3}, ps_nt)

@inferred n1(ps_nt)
@code_warntype n1(ps_nt)

@generated function (m::nphasetype1{T,K})(ps) where {T,K}
    
    # kv = :([:(MT.from_nt(im[$i], ps)) for i in 1:N])
    kv = [:(MT.from_nt(T.parameters[$i].parameters[1], ps)) for i in 1:length(T.types)]
    expr = :(Tuple([$(kv...)]))
    @show expr
    return expr
end

@inferred n1(ps_nt)


@generated function helper_make_nt2(m_tup::NTuple{N, K}, ps) where {N, K}
    
    # kv = :([:(MT.from_nt(im[$i], ps)) for i in 1:N])
    kv = [:(MT.from_nt(m_tup[$i], ps)) for i in 1:N]
    return kv
end

T_ = Tuple{SEO3, Ni2011}
@inferred n1(ps_nt)
@code_warntype n1(ps_nt)


mo1 = n1(ps_nt)
mutable struct nphasemodel{T <: nphasetype{nT, nK}, K, T1}
    m_vec::T
    mix::K
    ϕ::T1
end


n2 = two_phase_modelType(SEO3, Ni2011, mix_)

typeof(n2)
@inferred n2(ps_nt)

k = Type{SEO3}



using BenchmarkTools

@benchmark n1(ps_nt)
@benchmark n2(ps_nt)