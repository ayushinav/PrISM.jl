function get_Z(ρ::T1, h::T2, ω::T) where {T1, T2, T}
    broadcast!(exp10, ρ, ρ)
    k = sqrt(im * ω * μ / ρ[end])
    Z = ω * μ / k

    j = length(h)
    @inbounds while j >= 1
        k = sqrt(im * ω * μ / ρ[j])
        Z = ω * μ / k * coth(-im * k * h[j] + acoth(Z / (ω * μ / k)))
        j -= 1
    end

    broadcast!(log10, ρ, ρ)
    Z = conj(Z)
    return get_appres(Z, ω), get_phase(Z)
end

# the following are defined on scalars, there in-place variants don't make sense

"""
`get_phase(Z)`: returns the phase for impedance
"""
get_phase(Z) = 180 / π * atan(imag(Z) / real(Z));

"""
`get_appres(Z, ω)`: returns the ρₐ for impedance
"""
get_appres(Z, ω) = abs(Z)^2 / (eltype(ω))(μ) / ω;
