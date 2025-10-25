SubsurfaceCore.get_scales(::Type{<:MTResponse}, ::Val{:ρₐ}) = log10, log10
SubsurfaceCore.get_scales(::Type{<:MTResponse}, ::Val{:ϕ}) = log10, identity

SubsurfaceCore.get_labels(::Type{<:MTResponse}, ::Val{:ρₐ}) = "T(s)", "ρₐ (Ωm)"
SubsurfaceCore.get_labels(::Type{<:MTResponse}, ::Val{:ϕ}) = "T(s)", "ϕ (ᴼ)"

SubsurfaceCore.get_scales(::Type{<:MTModel}) = log10, log10
SubsurfaceCore.get_labels(::Type{<:MTModel}) = "ρ (Ωm)", "h(m)"