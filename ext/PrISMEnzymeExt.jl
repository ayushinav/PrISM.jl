module PrISMEnzymeExt

using PrISM, Enzyme
import .EnzymeRules: forward, reverse, augmented_primal
using .EnzymeRules
import PrISM: _dltar_c, get_c!

function EnzymeRules.forward(::FwdConfigWidth{N}, func::Const{typeof(get_c!)},
        ::Type{RT}, resp::Enzyme.Annotation, t::Enzyme.Annotation,
        m::Enzyme.Annotation{Tm}, mode::Enzyme.Annotation, dc::Enzyme.Annotation,
        c1::Enzyme.Annotation, c2::Enzyme.Annotation) where {RT, Tm, N}
    func.val(resp.val, t.val, m.val, mode.val, dc.val, c1.val, c2.val)

    ε = cbrt(eps(first(resp.val)))

    m1 = deepcopy(m.val)
    m2 = deepcopy(m.val)

    for i in eachindex(t.val) # this can be parallelized
        ω = 2π / t.val[i]
        fₓ = (_dltar_c(resp.val[i] + ε, ω, m.val) - _dltar_c(resp.val[i] - ε, ω, m.val)) /
             (2ε)
        for n in 1:N
            resp.dval[n][i] = 0
        end
        for k in fieldnames(Tm)
            for j in eachindex(getfield(m.val, k))
                for n in 1:N
                    k_dval = getfield(m.dval[n], k)
                    k_dval === getfield(m.val, k) && continue
                    getfield(m1, k)[j] += ε * k_dval[j]
                    getfield(m2, k)[j] -= ε * k_dval[j]
                    fₚⱼ = (_dltar_c(resp.val[i], ω, m1) - _dltar_c(resp.val[i], ω, m2)) /
                          (2ε)
                    getfield(m1, k)[j] -= ε * k_dval[j]
                    getfield(m2, k)[j] += ε * k_dval[j]
                    resp.dval[n][i] += (-fₚⱼ / fₓ)
                end
            end
        end
    end
    return nothing
end

function EnzymeRules.augmented_primal(config::RevConfigWidth, func::Const{typeof(get_c!)},
        ::Type{RT}, resp::Enzyme.Annotation, t::Enzyme.Annotation,
        m::Enzyme.Annotation{Tm}, mode::Enzyme.Annotation, dc::Enzyme.Annotation,
        c1::Enzyme.Annotation, c2::Enzyme.Annotation) where {RT, Tm}
    func.val(resp.val, t.val, m.val, mode.val, dc.val, c1.val, c2.val)
    primal = nothing
    tape = overwritten(config)[4] ? deepcopy(m.val) : nothing
    return AugmentedReturn(primal, nothing, tape)
end

function EnzymeRules.reverse(config::RevConfigWidth{1}, func::Const{typeof(get_c!)},
        ::Type{RT}, tape, resp::Enzyme.Annotation, t::Enzyme.Annotation,
        m::Enzyme.Annotation{Tm}, mode::Enzyme.Annotation, dc::Enzyme.Annotation,
        c1::Enzyme.Annotation, c2::Enzyme.Annotation) where {RT, Tm}
    m_val = overwritten(config)[4] ? tape : m.val

    ε = cbrt(eps(first(resp.val)))

    m1 = deepcopy(m_val)
    m2 = deepcopy(m_val)

    for i in eachindex(t.val) # this can be parallelized
        ω = 2π / t.val[i]
        fₓ = (_dltar_c(resp.val[i] + ε, ω, m_val) - _dltar_c(resp.val[i] - ε, ω, m_val)) /
             (2ε)
        for k in fieldnames(Tm)
            k_val = getfield(m.val, k)
            k_dval = getfield(m.dval, k)
            k_dval === k_val && continue
            for j in eachindex(getfield(m_val, k))
                getfield(m1, k)[j] += ε
                getfield(m2, k)[j] -= ε
                fₚⱼ = (_dltar_c(resp.val[i], ω, m1) - _dltar_c(resp.val[i], ω, m2)) / (2ε)
                k_dval[j] += resp.dval[i] * (-fₚⱼ / fₓ)
                getfield(m1, k)[j] -= ε
                getfield(m2, k)[j] += ε
            end
        end
        resp.dval[i] = 0
        # end
    end

    d_dc = typeof(dc) <: Active ? (zero(dc.val)) : nothing
    d_c1 = typeof(c1) <: Active ? (zero(c1.val)) : nothing
    d_c2 = typeof(c2) <: Active ? (zero(c2.val)) : nothing

    return (nothing, nothing, nothing, nothing, d_dc, d_c1, d_c2)
end

function EnzymeRules.reverse(config::RevConfigWidth{N}, func::Const{typeof(get_c!)},
        ::Type{RT}, tape, resp::Enzyme.Annotation, t::Enzyme.Annotation,
        m::Enzyme.Annotation{Tm}, mode::Enzyme.Annotation, dc::Enzyme.Annotation,
        c1::Enzyme.Annotation, c2::Enzyme.Annotation) where {RT, Tm, N}
    m_val = overwritten(config)[4] ? tape : m.val

    ε = cbrt(eps(first(resp.val)))

    m1 = deepcopy(m_val)
    m2 = deepcopy(m_val)

    for i in eachindex(t.val) # this can be parallelized
        ω = 2π / t.val[i]
        fₓ = (_dltar_c(resp.val[i] + ε, ω, m_val) - _dltar_c(resp.val[i] - ε, ω, m_val)) /
             (2ε)

        for k in fieldnames(Tm)
            for j in eachindex(getfield(m_val, k))
                getfield(m1, k)[j] += ε
                getfield(m2, k)[j] -= ε
                fₚⱼ = (_dltar_c(resp.val[i], ω, m1) - _dltar_c(resp.val[i], ω, m2)) / (2ε)
                getfield(m1, k)[j] -= ε
                getfield(m2, k)[j] += ε
                for n in 1:N
                    k_dval = getfield(m.dval[n], k)
                    k_dval === getfield(m.val, k) && continue
                    k_dval[j] += resp.dval[n][i] * (-fₚⱼ / fₓ)
                end
            end
        end
        for n in 1:N
            resp.dval[n][i] = 0
        end
    end

    d_dc = typeof(dc) <: Active ? ntuple(_ -> zero(dc.val), Val(N)) : nothing
    d_c1 = typeof(c1) <: Active ? ntuple(_ -> zero(c1.val), Val(N)) : nothing
    d_c2 = typeof(c2) <: Active ? ntuple(_ -> zero(c2.val), Val(N)) : nothing

    return (nothing, nothing, nothing, nothing, d_dc, d_c1, d_c2)
end

end
