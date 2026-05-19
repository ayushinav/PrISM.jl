module PrISMEnzymeExt

using PrISM, Enzyme
import .EnzymeRules: forward, reverse, augmented_primal
using .EnzymeRules
import PrISM: find_c, get_c!

function _ift_find_c_enzyme(f_func, c, ω_val, m_val, m_dval)
    # @show "rev"

    # c1 = one(c)
    # f1(c_) = f_func(c_, ω_val, m_val)
    # # fₓ = first(Enzyme.autodiff(Enzyme.Forward, Const(f1),
    # #                            Enzyme.Duplicated(c, c1)))
    # fₓ = DifferentiationInterface.derivative(f1, AutoEnzyme(; mode = Reverse), c)
    # # fₚ = first(Enzyme.autodiff(Enzyme.Forward, m_ -> f_func(c, ω_val, m_),
    # #                            Enzyme.Duplicated(m_val, m_dval)))
    # # return -fₚ / fₓ
    # return fₓ

    ε = sqrt(eps(c))
    f1 = f_func(c + ε, ω_val, m_val)
    f2 = f_func(c - ε, ω_val, m_val)
    fₓ = (f1 - f2) / (2ε)
    m_pert1 = deepcopy(m_val)
    m_pert2 = deepcopy(m_val)

    m_pert1.m .= m_val.m .+ ε .* m_dval.m
    m_pert1.h .= m_val.h .+ ε .* m_dval.h
    m_pert1.ρ .= m_val.ρ .+ ε .* m_dval.ρ
    m_pert1.vp .= m_val.vp .+ ε .* m_dval.vp

    m_pert2.m .= m_val.m .- ε .* m_dval.m
    m_pert2.h .= m_val.h .- ε .* m_dval.h
    m_pert2.ρ .= m_val.ρ .- ε .* m_dval.ρ
    m_pert2.vp .= m_val.vp .- ε .* m_dval.vp

    fₚ = (f_func(c, ω_val, m_pert1) - f_func(c, ω_val, m_pert2)) / (2ε)
    -fₚ / fₓ
end

# function _ift_find_c_enzyme2(f_func, c, ω_val, m)
#     @show typeof(c), typeof(ω_val), typeof(m)

#     c1 = one(c)
#     f1(c_) = f_func(c_, ω_val, m.val)
#     fₓ = first(Enzyme.autodiff(Enzyme.Forward, Const(f1), Enzyme.Duplicated(c, c1)))
#     fₚ = first(Enzyme.autodiff(Enzyme.Forward, m_ -> f_func(c, ω_val, m_), m))
#     return -fₚ / fₓ
#     # return fₓ
#     # return 1.
# end

# RWModel

function EnzymeRules.forward(::FwdConfigWidth, ::Const{typeof(find_c)},
        ::Type{<:Duplicated}, f::Const, c1::T1, c2::T2, ω::Enzyme.Const, dc::Enzyme.Const,
        m::Duplicated{<:RWModel}) where {T1 <: Enzyme.Annotation, T2 <: Enzyme.Annotation}
    println("Using custom rule! : Duplicated : unbatched")

    c = find_c(f.val, c1.val, c2.val, ω.val, dc.val, m.val)
    d_c = _ift_find_c_enzyme(f.val, c, ω.val, m.val, m.dval)
    # dc = _ift_find_c_enzyme2(f.val, c, ω.val, m)

    return Duplicated(c, d_c)
end

function EnzymeRules.forward(::FwdConfigWidth, ::Const{typeof(find_c)},
        ::Type{<:DuplicatedNoNeed}, f::Const, c1::T1, c2::T2, ω::Enzyme.Const, dc::Enzyme.Const,
        m::Duplicated{<:RWModel}) where {T1 <: Enzyme.Annotation, T2 <: Enzyme.Annotation}
    println("Using custom rule! : DuplicatedNoNeed : unbatched")

    c = find_c(f.val, c1.val, c2.val, ω.val, dc.val, m.val)
    d_c = _ift_find_c_enzyme(f.val, c, ω.val, m.val, m.dval)
    # dc = _ift_find_c_enzyme2(f.val, c, ω.val, m.val)

    return d_c
end

# batch: fₓ computed once, fₚ per direction
function _ift_find_c_enzyme_batch(f_func, c, ω_val, m_val, m_dvals::NTuple{N}) where {N}
    ε = sqrt(eps(c))
    f1 = f_func(c + ε, ω_val, m_val)
    f2 = f_func(c - ε, ω_val, m_val)
    fₓ = (f1 - f2) / (2ε)
    m_pert1 = deepcopy(m_val)
    m_pert2 = deepcopy(m_val)
    return ntuple(Val(N)) do i
        # fₚ = first(Enzyme.autodiff(Enzyme.Forward, m_ -> f_func(c, ω_val, m_),
        #                            Enzyme.Duplicated(m_val, m_dvals[i])))
        # -fₚ / fₓ
        # fp1 = f_fun()

        md = m_dvals[i]
        m_pert1.m .= m_val.m .+ ε .* md.m
        m_pert1.h .= m_val.h .+ ε .* md.h
        m_pert1.ρ .= m_val.ρ .+ ε .* md.ρ
        m_pert1.vp .= m_val.vp .+ ε .* md.vp

        m_pert2.m .= m_val.m .- ε .* md.m
        m_pert2.h .= m_val.h .- ε .* md.h
        m_pert2.ρ .= m_val.ρ .- ε .* md.ρ
        m_pert2.vp .= m_val.vp .- ε .* md.vp

        fₚ = (f_func(c, ω_val, m_pert1) - f_func(c, ω_val, m_pert2)) / (2ε)
        -fₚ / fₓ
    end
end

function EnzymeRules.forward(::FwdConfigWidth,
        ::Const{typeof(find_c)},
        ::Type{<:BatchDuplicated},
        f::Const,
        c1::T1,
        c2::T2,
        ω::Enzyme.Const, dc::Enzyme.Const,
        m::BatchDuplicated{<:RWModel, N}) where {
        T1 <: Enzyme.Annotation, T2 <: Enzyme.Annotation, N}
    println("Using custom rule! : BatchDuplicated : batched")
    c = find_c(f.val, c1.val, c2.val, ω.val, dc.val, m.val)
    dcs = _ift_find_c_enzyme_batch(f.val, c, ω.val, m.val, m.dval)
    return BatchDuplicated(c, dcs)
end

function EnzymeRules.forward(::FwdConfigWidth,
        ::Const{typeof(find_c)},
        ::Type{<:BatchDuplicatedNoNeed},
        f::Const,
        c1::T1,
        c2::T2,
        ω::Enzyme.Const, dc::Enzyme.Const,
        m::BatchDuplicated{<:RWModel, N}) where {
        T1 <: Enzyme.Annotation, T2 <: Enzyme.Annotation, N}
    println("Using custom rule! : BatchDuplicatedNoNeed : batched")
    c = find_c(f.val, c1.val, c2.val, ω.val, dc.val, m.val)
    dcs = _ift_find_c_enzyme_batch(f.val, c, ω.val, m.val, m.dval)
    return dcs
end



















function _ift_find_c_enzyme(f_func, c, ω_val, m_val::LWModel, m_dval)
    # @show typeof(c), typeof(ω_val), typeof(m_val), typeof(m_dval)
    # @show "rev"

    # c1 = one(c)
    # f1(c_) = f_func(c_, ω_val, m_val)
    # # fₓ = first(Enzyme.autodiff(Enzyme.Forward, Const(f1),
    # #                            Enzyme.Duplicated(c, c1)))
    # fₓ = DifferentiationInterface.derivative(f1, AutoEnzyme(; mode = Reverse), c)
    # # fₚ = first(Enzyme.autodiff(Enzyme.Forward, m_ -> f_func(c, ω_val, m_),
    # #                            Enzyme.Duplicated(m_val, m_dval)))
    # # return -fₚ / fₓ
    # return fₓ

    ε = sqrt(eps(c))
    f1 = f_func(c + ε, ω_val, m_val)
    f2 = f_func(c - ε, ω_val, m_val)
    fₓ = (f1 - f2) / (2ε)
    m_pert1 = deepcopy(m_val)
    m_pert2 = deepcopy(m_val)

    m_pert1.m .= m_val.m .+ ε .* m_dval.m
    m_pert1.h .= m_val.h .+ ε .* m_dval.h
    m_pert1.ρ .= m_val.ρ .+ ε .* m_dval.ρ

    m_pert2.m .= m_val.m .- ε .* m_dval.m
    m_pert2.h .= m_val.h .- ε .* m_dval.h
    m_pert2.ρ .= m_val.ρ .- ε .* m_dval.ρ

    fₚ = (f_func(c, ω_val, m_pert1) - f_func(c, ω_val, m_pert2)) / (2ε)
    -fₚ / fₓ
    # return 1.
end

#LWModel

function EnzymeRules.forward(::FwdConfigWidth, ::Const{typeof(find_c)},
        ::Type{<:Duplicated}, f::Const, c1::T1, c2::T2, ω::Enzyme.Const, dc::Enzyme.Const,
        m::Duplicated{<:LWModel}) where {T1 <: Enzyme.Annotation, T2 <: Enzyme.Annotation}
    # println("Using custom rule! : Duplicated")

    c = find_c(f.val, c1.val, c2.val, ω.val, dc.val, m.val)
    d_c = _ift_find_c_enzyme(f.val, c, ω.val, m.val, m.dval)
    # dc = _ift_find_c_enzyme2(f.val, c, ω.val, m)

    return Duplicated(c, d_c)
end

function EnzymeRules.forward(::FwdConfigWidth, ::Const{typeof(find_c)},
        ::Type{<:DuplicatedNoNeed}, f::Const, c1::T1, c2::T2, ω::Enzyme.Const, dc::Enzyme.Const,
        m::Duplicated{<:LWModel}) where {T1 <: Enzyme.Annotation, T2 <: Enzyme.Annotation}
    # println("Using custom rule! : DuplicatedNoNeed")

    c = find_c(f.val, c1.val, c2.val, ω.val, dc.val, m.val)
    d_c = _ift_find_c_enzyme(f.val, c, ω.val, m.val, m.dval)
    # dc = _ift_find_c_enzyme2(f.val, c, ω.val, m.val)

    return d_c
end

# batch: fₓ computed once, fₚ per direction
function _ift_find_c_enzyme_batch(
        f_func, c, ω_val, m_val::LWModel, m_dvals::NTuple{N}) where {N}
    ε = sqrt(eps(c))
    f1 = f_func(c + ε, ω_val, m_val)
    f2 = f_func(c - ε, ω_val, m_val)
    fₓ = (f1 - f2) / (2ε)
    m_pert1 = deepcopy(m_val)
    m_pert2 = deepcopy(m_val)
    return ntuple(Val(N)) do i
        # fₚ = first(Enzyme.autodiff(Enzyme.Forward, m_ -> f_func(c, ω_val, m_),
        #                            Enzyme.Duplicated(m_val, m_dvals[i])))
        # -fₚ / fₓ
        # fp1 = f_fun()

        md = m_dvals[i]
        m_pert1.m .= m_val.m .+ ε .* md.m
        m_pert1.h .= m_val.h .+ ε .* md.h
        m_pert1.ρ .= m_val.ρ .+ ε .* md.ρ

        m_pert2.m .= m_val.m .- ε .* md.m
        m_pert2.h .= m_val.h .- ε .* md.h
        m_pert2.ρ .= m_val.ρ .- ε .* md.ρ

        fₚ = (f_func(c, ω_val, m_pert1) - f_func(c, ω_val, m_pert2)) / (2ε)
        -fₚ / fₓ
    end
end

function EnzymeRules.forward(::FwdConfigWidth,
        ::Const{typeof(find_c)},
        ::Type{<:BatchDuplicated},
        f::Const,
        c1::T1,
        c2::T2,
        ω::Enzyme.Const, dc::Enzyme.Const,
        m::BatchDuplicated{<:LWModel, N}) where {
        T1 <: Enzyme.Annotation, T2 <: Enzyme.Annotation, N}
    # println("Using custom rule! : BatchDuplicated 2")
    c = find_c(f.val, c1.val, c2.val, ω.val, dc.val, m.val)
    dcs = _ift_find_c_enzyme_batch(f.val, c, ω.val, m.val, m.dval)
    return BatchDuplicated(c, dcs)
end

function EnzymeRules.forward(::FwdConfigWidth,
        ::Const{typeof(find_c)},
        ::Type{<:BatchDuplicatedNoNeed},
        f::Const,
        c1::T1,
        c2::T2,
        ω::Enzyme.Const, dc::Enzyme.Const,
        m::BatchDuplicated{<:LWModel, N}) where {
        T1 <: Enzyme.Annotation, T2 <: Enzyme.Annotation, N}
    c = find_c(f.val, c1.val, c2.val, ω.val, dc.val, m.val)
    dcs = _ift_find_c_enzyme_batch(f.val, c, ω.val, m.val, m.dval)
    return dcs
end

function EnzymeRules.augmented_primal(config::RevConfigWidth, func::Const{typeof(get_c!)},
        ::Type{RT}, resp::Enzyme.Annotation, t::Enzyme.Annotation,
        m::Enzyme.Annotation{Tm}, mode::Enzyme.Annotation, dc::Enzyme.Annotation,
        c1::Enzyme.Annotation, c2::Enzyme.Annotation) where {RT, Tm}
    println("Agumented primal")

    func.val(resp.val, t.val, m.val, mode.val, dc.val, c1.val, c2.val)
    primal = nothing
    tape = deepcopy(m.val) #overwritten(config)[4] ? deepcopy(m.val) : nothing
    println("h[1] in augmented_primal: ", m.val.h[1])  # should be ~0.1 (km), not ~100 (m)
    return AugmentedReturn(primal, nothing, tape)
end

function EnzymeRules.reverse(config::RevConfigWidth{1}, func::Const{typeof(get_c!)},
        ::Type{RT}, tape, resp::Enzyme.Annotation, t::Enzyme.Annotation,
        m::Enzyme.Annotation{Tm}, mode::Enzyme.Annotation, dc::Enzyme.Annotation,
        c1::Enzyme.Annotation, c2::Enzyme.Annotation) where {RT, Tm}

    # @show tape.h
    # @show m.val.h
    # @show typeof(m)
    println("only m")
    m_val = overwritten(config)[4] ? tape : m.val

    @show m_val.h
    ε = sqrt(eps(first(resp.val)))

    m1 = deepcopy(m_val)
    m2 = deepcopy(m_val)

    for i in eachindex(t.val) # this can be parallelized
        ω = 2π / t.val[i]
        fₓ = (_dltar_c(resp.val[i] + ε, ω, m_val) - _dltar_c(resp.val[i] - ε, ω, m_val)) /
             (2ε)

        # for im in 1:(mode + 1)

        # for (field, shadow_field) in ((:m, :m), (:h, :h), (:ρ, :ρ), (:vp, :vp))
        #     sfld = getfield(m.dval, shadow_field)
        #     for j in eachindex(getfield(mval, field))
        #         getfield(m1, field)[j] += ε
        #         getfield(m2, field)[j] -= ε
        #         fₚⱼ = (f.val(c, ω.val, m1) - f.val(c, ω.val, m2)) / (2ε)
        #         sfld[j] += dret.val * (-fₚⱼ / fₓ)
        #         getfield(m1, field)[j] -= ε
        #         getfield(m2, field)[j] += ε
        #     end
        # end

        for k in fieldnames(Tm)
            # k_dval = getfield(m.dval, k)
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
        # end
    end

    d_dc = typeof(c2) <: Active ? (zero(dc.val)) : nothing
    d_c1 = typeof(c1) <: Active ? (zero(c1.val)) : nothing
    d_c2 = typeof(c2) <: Active ? (zero(c2.val)) : nothing

    return (nothing, nothing, nothing, nothing, d_dc, d_c1, d_c2)
end

function EnzymeRules.reverse(config::RevConfigWidth{N}, func::Const{typeof(get_c!)},
        ::Type{RT}, tape, resp::Enzyme.Annotation, t::Enzyme.Annotation,
        m::Enzyme.Annotation{Tm}, mode::Enzyme.Annotation, dc::Enzyme.Annotation,
        c1::Enzyme.Annotation, c2::Enzyme.Annotation) where {RT, Tm, N}

    # @show tape.h
    # @show m.val.h
    # @show typeof(m)
    println("only m")
    m_val = overwritten(config)[4] ? tape : m.val

    @show m_val.h
    @show resp.val
    ε = sqrt(eps(first(resp.val)))

    m1 = deepcopy(m_val)
    m2 = deepcopy(m_val)

    for i in eachindex(t.val) # this can be parallelized
        ω = 2π / t.val[i]
        fₓ = (_dltar_c(resp.val[i] + ε, ω, m_val) - _dltar_c(resp.val[i] - ε, ω, m_val)) /
             (2ε)

        # for im in 1:(mode + 1)

        # for (field, shadow_field) in ((:m, :m), (:h, :h), (:ρ, :ρ), (:vp, :vp))
        #     sfld = getfield(m.dval, shadow_field)
        #     for j in eachindex(getfield(mval, field))
        #         getfield(m1, field)[j] += ε
        #         getfield(m2, field)[j] -= ε
        #         fₚⱼ = (f.val(c, ω.val, m1) - f.val(c, ω.val, m2)) / (2ε)
        #         sfld[j] += dret.val * (-fₚⱼ / fₓ)
        #         getfield(m1, field)[j] -= ε
        #         getfield(m2, field)[j] += ε
        #     end
        # end

        for k in fieldnames(Tm)
            # k_dval = getfield(m.dval, k)
            # k_val = getfield(m.val, k)
            # k_dval = getfield(m[n].dval, k)
            # k_dval === k_val && continue
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
    end

    d_dc = typeof(dc) <: Active ? ntuple(_ -> zero(dc.val), Val(N)) : nothing
    d_c1 = typeof(c1) <: Active ? ntuple(_ -> zero(c1.val), Val(N)) : nothing
    d_c2 = typeof(c2) <: Active ? ntuple(_ -> zero(c2.val), Val(N)) : nothing


    return (nothing, nothing, nothing, nothing, d_dc, d_c1, d_c2)
end

end