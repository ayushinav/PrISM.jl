function Base.show(io::IO, m::model) where {model <:
                                            MTModel{<:AbstractVector, <:AbstractVector}}
    println("1D ", typeof(m).name.name, " : ")
    println("Layer \t h \t log(ρ)")
    println("____________________________")
    for i in eachindex(m.h)
        m_ = round(m.m[i]; digits=3)
        h_ = round(m.h[i]; digits=3)
        println("$i \t $h_ \t $m_")
    end
    m_ = round(m.m[end]; digits=3)
    println("$(length(m.m)) \t ∞ \t \t $m_")
end

function Base.show(io::IO,
        m::model) where {model <: RWModel{
        <:AbstractVector, <:AbstractVector, <:AbstractVector, <:AbstractVector}}
    println("1D ", typeof(m).name.name, " : ")
    println("Layer \t h(km) \t Vₛ (km/s) \t ρ (kg⋅m⁻³) \t Vₚ (km/s)")
    println("_________________________________________________")
    for i in eachindex(m.h)
        m_ = round(m.m[i]; digits=3)
        h_ = round(m.h[i]; digits=3)
        d_ = round(m.ρ[i]; digits=3)
        vp_ = round(m.vp[i]; digits=3)
        println("$i \t $h_ \t $m_ \t $d_ \t $vp_")
    end
    m_ = round(m.m[end]; digits=3)
    d_ = round(m.ρ[end]; digits=3)
    vp_ = round(m.vp[end]; digits=3)
    println("$(length(m.m)) \t ∞ \t \t $m_ \t $d_ \t $vp_")
end

function Base.show(io::IO,
        m::model) where {model <:
                         LWModel{<:AbstractVector, <:AbstractVector, <:AbstractVector}}
    println("1D ", typeof(m).name.name, " : ")
    println("Layer  h(km) \t Vₛ (km/s) \t  ρ (kg⋅m⁻³)")
    println("_________________________________________________")
    for i in eachindex(m.h)
        m_ = round(m.m[i]; digits=3)
        h_ = round(m.h[i]; digits=3)
        d_ = round(m.ρ[i]; digits=3)
        println("$i \t $h_ \t $m_ \t $d_")
    end
    m_ = round(m.m[end]; digits=3)
    d_ = round(m.ρ[end]; digits=3)
    println("$(length(m.m)) \t ∞ \t \t $m_ \t $d_")
end

function Base.show(io::IO, m::model) where {model <:
                                            DCModel{<:AbstractVector, <:AbstractVector}}
    println("1D ", typeof(m).name.name, " : ")
    println("Layer \t h \t log(ρ)")
    println("____________________________")
    for i in eachindex(m.h)
        m_ = round(m.m[i]; digits=3)
        h_ = round(m.h[i]; digits=3)
        println("$i \t $h_ \t $m_")
    end
    m_ = round(m.m[end]; digits=3)
    println("$(length(m.m)) \t ∞ \t \t $m_")
end