SubsurfaceCore.get_scales(::Type{<:MTResponse}, ::Val{:ρₐ}) = log10, log10
SubsurfaceCore.get_scales(::Type{<:MTResponse}, ::Val{:ϕ}) = log10, identity

SubsurfaceCore.get_labels(::Type{<:MTResponse}, ::Val{:ρₐ}) = "T(s)", "ρₐ (Ωm)"
SubsurfaceCore.get_labels(::Type{<:MTResponse}, ::Val{:ϕ}) = "T(s)", "ϕ (ᴼ)"

SubsurfaceCore.get_scales(::Type{<:MTModel}) = log10, log10
SubsurfaceCore.get_labels(::Type{<:MTModel}) = "ρ (Ωm)", "h(m)"

SubsurfaceCore.get_scales(::Type{<:SurfaceWaveResponse}, ::Val{:c}) = log10, identity

function SubsurfaceCore.get_labels(::Type{<:SurfaceWaveResponse}, ::Val{:c})
    "T(s)", "Velocity (km/s)"
end

function SubsurfaceCore.get_scales(::Type{T}) where {T <: Union{LWModel, RWModel}}
    identity, identity
end
function SubsurfaceCore.get_labels(::Type{T}) where {T <: Union{LWModel, RWModel}}
    "vₛ (km/s)", "h(km)"
end
