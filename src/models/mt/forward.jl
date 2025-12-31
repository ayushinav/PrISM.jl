const global μ = 4π * 1.0f-7; # Float32 will promote to Float64 without a problem

"""
`get_Z(ρ,h,ω)`:

returns a tuple of ρₐ and ϕ, given arrays of resistivity `ρ` and thickness `h` for the angular frequenciy `ω`.
"""

const default_mt_tf_fns = (ρₐ=no_tf, ϕ=no_tf)

# dispatch on forward for 1d model
"""
`forward(m::model, ω::Vector{T}) where T <: Union{Float32, Float64}`:

returns a  `response` for the given model `m` at the frequencies  `ω`
"""
function SubsurfaceCore.forward(m::Tm, ω::T3, params=default_params_mt) where {Tm <: MTModel, T3}
    if !(length(m.h) == length(m.m) - 1)
        error("number of model layers should be 1 less than the number of model parameters")
    end
    n = length(ω)
    ρₐ = zeros(eltype(m.m), n)
    ϕ = zeros(eltype(m.m), n)
    i = 1
    @inbounds while i <= n
        ρₐ[i], ϕ[i] = get_Z(m.m, m.h, ω[i])
        i += 1
    end
    MTResponse(ρₐ, ϕ)
end

# dispatch on forward! for 1d model

"""
`forward!(r, m, ω)

    updates response `r` type for the given model `m` at the frequencies  `ω`
"""
function forward!(r::Tr, m::Tm, ω::T3, params=default_params_mt) where {Tr <: MTResponse, Tm <: MTModel, T3}
    if !(length(m.h) == length(m.m) - 1)
        error("number of model layers should be 1 less than the number of model parameters")
    end
    n = length(ω)
    i = 1
    @inbounds while i <= n
        r.ρₐ[i], r.ϕ[i] = get_Z(m.m, m.h, ω[i])
        i += 1
    end
    nothing
end
