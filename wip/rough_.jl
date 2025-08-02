function from_nt(m::Type{T}, nt::NamedTuple) where {T <: two_phase_model}
    ϕ = nt.ϕ
    m1 = T.parameters[2]
    m2 = T.parameters[3]
    mix = T.parameters[4]

    model1 = MT.from_nt(m1, nt)
    model2 = MT.from_nt(m2, nt)

    return two_phase_model(ϕ, model1, model2, mix())
end