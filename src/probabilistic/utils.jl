import SubsurfaceCore: sample_type

for m in subtypes(AbstractGeophyModel)
    mstring = string(m)
    dstring = mstring * "Distribution"
    eval(Meta.parse("SubsurfaceCore.sample_type(d::$dstring) = $mstring"))
    eval(Meta.parse("SubsurfaceCore.sample_type(::Type{T}) where {T <:$dstring} = $mstring"))
end
