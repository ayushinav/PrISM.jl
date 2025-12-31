mutable struct DCResponse{T1} <: AbstractGeophyResponse
    ρₐ::T1
end

#=
More often than not, T1 and T2 will be same, that is typeof(ρₐ) and typeof(ϕ) will be same but this might become different when using a rock physics model in front of this. 
For now, we go with the same design as for `DCModel`.
=#

"""
create a model type for a given resistivity distribution
that can be used to calculate forward response for 1d DC

## Usage

```jldoctest
m = [4.0, 2.0, 3.0]
h = [1000.0, 1000.0]
model = DCModel(m, h)
print(model)

# output

1D DCModel : 
Layer    h       log(ρ)
____________________________
1        1000.0          4.0
2        1000.0          2.0
3        ∞               3.0
```
"""
mutable struct DCModel{T1 <: AbstractArray{<:Any}, T2 <: AbstractArray{<:Any}} <:
               AbstractGeophyModel{T1, T2}
    m::T1
    h::T2
end

# explain in blog the reasoning behind this! This covers all 1D, 2D, 3D models for DC.
# make pretty tables to print these models

SubsurfaceCore.default_params(::Type{T}) where {T <: DCModel} = (; hankel_filter=Key101)
const default_params_DC = default_params(DCModel)
