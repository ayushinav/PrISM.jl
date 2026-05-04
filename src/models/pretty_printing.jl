function Base.show(io::IO, m::model) where {model <:
                                            MTModel{<:AbstractVector, <:AbstractVector}}
    
    mat_ = hcat(
        eachindex(m.m), 
        round.(m.m; digits = 3), 
        [round.(m.h; digits = 3)..., "∞"], 
        [round.(cumsum(m.h); digits = 3)..., "∞"])
    
    header_ = (
        ["Layer",  "log-ρ", "h", "z"],
        ["", "[Ωm]", "[m]", "[m]"]
    )
    pretty_table(io, mat_; 
                header = header_,
                formatters    = ft_printf("%5.2f", 2:4),
                header_crayon = crayon"blue bold",
                title = "1D MT Model")
end

function Base.show(io::IO,
        m::model) where {model <: RWModel{
        <:AbstractVector, <:AbstractVector, <:AbstractVector, <:AbstractVector}}
    
    mat_ = hcat(
        eachindex(m.m), 
        round.(m.m; digits = 3), 
        [round.(m.h; digits = 3)..., "∞"], 
        [round.(cumsum(m.h); digits = 3)..., "∞"],
        round.(m.ρ; digits = 3), 
        round.(m.vp; digits = 3))
    
    header_ = (
        ["Layer",  "vₛ", "h", "z", "ρ", "vₚ"],
        ["", "[km/s]", "[m]", "[m]", "[kg⋅m⁻³]", "[km/s]"]
    )
    pretty_table(io, mat_; 
                header = header_,
                formatters    = ft_printf("%5.2f", 2:4),
                header_crayon = crayon"blue bold",
                title = "1D Rayleigh Wave Model")
end

function Base.show(io::IO,
        m::model) where {model <:
                         LWModel{<:AbstractVector, <:AbstractVector, <:AbstractVector}}
    mat_ = hcat(
        eachindex(m.m), 
        round.(m.m; digits = 3), 
        [round.(m.h; digits = 3)..., "∞"], 
        [round.(cumsum(m.h); digits = 3)..., "∞"],
        round.(m.ρ; digits = 3))
    
    header_ = (
        ["Layer",  "vₛ", "h", "z", "ρ"],
        ["", "[km/s]", "[m]", "[m]", "[kg⋅m⁻³]"]
    )
    pretty_table(io, mat_; 
                header = header_,
                formatters    = ft_printf("%5.2f", 2:4),
                header_crayon = crayon"blue bold",
                title = "1D Love Wave Model")
end

function Base.show(io::IO, m::model) where {model <:
                                            DCModel{<:AbstractVector, <:AbstractVector}}
    
    mat_ = hcat(
        eachindex(m.m), 
        round.(m.m; digits = 3), 
        [round.(m.h; digits = 3)..., "∞"], 
        [round.(cumsum(m.h); digits = 3)..., "∞"])
    
    header_ = (
        ["Layer",  "log-ρ", "h", "z"],
        ["", "[Ωm]", "[m]", "[m]"]
    )
    pretty_table(io, mat_; 
                header = header_,
                formatters    = ft_printf("%5.2f", 2:4),
                header_crayon = crayon"blue bold",
                title = "1D DC Resistivity Model")
end
