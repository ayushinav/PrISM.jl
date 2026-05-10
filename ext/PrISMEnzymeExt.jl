module PrISMEnzymeExt

using PrISM, Enzyme
import .EnzymeRules: forward, reverse, augmented_primal
using .EnzymeRules
import PrISM: find_c_

function _ift_find_c_enzyme(f_func, c, ω_val, m_val, m_dval)
    fₓ = first(Enzyme.autodiff(Enzyme.Forward, c_ -> f_func(c_, ω_val, m_val),
                               Enzyme.Duplicated(c, one(c))))
    fₚ = first(Enzyme.autodiff(Enzyme.Forward, m_ -> f_func(c, ω_val, m_),
                               Enzyme.Duplicated(m_val, m_dval)))
    return -fₚ / fₓ
end

# RWModel

function EnzymeRules.forward(
        ::FwdConfigWidth, ::Const{typeof(find_c_)}, 
        ::Type{<:Duplicated}, f::Const,
        c1::T1, c2::T2, ω::Enzyme.Const, m::Duplicated{<:RWModel}) where {T1 <: Enzyme.Annotation, T2 <: Enzyme.Annotation}
    println("Using custom rule!")

    c  = find_c_(f.val, c1.val, c2.val, ω.val, m.val)
    dc = _ift_find_c_enzyme(f.val, c, ω.val, m.val, m.dval)

    return Duplicated(c, dc)
end

function EnzymeRules.forward(
        ::FwdConfigWidth, ::Const{typeof(find_c_)}, 
        ::Type{<:DuplicatedNoNeed}, f::Const,
        c1::T1, c2::T2, ω::Enzyme.Const, m::Duplicated{<:RWModel}) where {T1 <: Enzyme.Annotation, T2 <: Enzyme.Annotation}
    println("Using custom rule!")

    c = find_c_(f.val, c1.val, c2.val, ω.val, m.val)
    dc = _ift_find_c_enzyme(f.val, c, ω.val, m.val, m.dval)

    return dc
end


# batch: fₓ computed once, fₚ per direction
function _ift_find_c_enzyme(f_func, c, ω_val, m_val, m_dvals::NTuple{N}) where {N}
    fₓ = only(Enzyme.autodiff(Enzyme.Forward, c_ -> f_func(c_, ω_val, m_val),
                               Enzyme.Duplicated(c, one(c))))
    return ntuple(N) do i
        fₚ = only(Enzyme.autodiff(Enzyme.Forward, m_ -> f_func(c, ω_val, m_),
                                   Enzyme.Duplicated(m_val, m_dvals[i])))
        # -fₚ / fₓ
        fₓ
    end
end

function EnzymeRules.forward(
        ::FwdConfigWidth, ::Const{typeof(find_c_)},
        ::Type{<:BatchDuplicated}, f::Const,
        c1::T1, c2::T2, ω::Enzyme.Const,
        m::BatchDuplicated{<:RWModel, N}) where {T1 <: Enzyme.Annotation, T2 <: Enzyme.Annotation, N}
    println("Using custom rule!")
    c   = find_c_(f.val, c1.val, c2.val, ω.val, m.val)
    dcs = _ift_find_c_enzyme(f.val, c, ω.val, m.val, m.dval)
    return BatchDuplicated(c, dcs)
end


function EnzymeRules.forward(
        ::FwdConfigWidth, ::Const{typeof(find_c_)},
        ::Type{<:BatchDuplicatedNoNeed}, f::Const,
        c1::T1, c2::T2, ω::Enzyme.Const,
        m::BatchDuplicated{<:RWModel, N}) where {T1 <: Enzyme.Annotation, T2 <: Enzyme.Annotation, N}
    println("Using custom rule!")
    c   = find_c_(f.val, c1.val, c2.val, ω.val, m.val)
    dcs = _ift_find_c_enzyme(f.val, c, ω.val, m.val, m.dval)
    return dcs
end

end