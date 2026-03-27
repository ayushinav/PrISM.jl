mutable struct SurfaceWaveResponse{T} <: AbstractGeophyResponse
    c::T
end
"""
create a model type for a given resistivity distribution
that can be used to calculate forward response for 1d MT

## Usage

```jldoctest
m = [4e3, 3e3, 3.5e3]
h = [1000.0, 1000.0]
d = [3e3, 3e3, 3e3]
model = LWModel(m, h, d)

# output

1D Love Wave Model
┌───────┬─────────┬─────────┬─────────┬──────────┐
│ Layer │      vₛ │       h │       z │        ρ │
│       │  [km/s] │     [m] │     [m] │ [kg⋅m⁻³] │
├───────┼─────────┼─────────┼─────────┼──────────┤
│     1 │ 4000.00 │ 1000.00 │ 1000.00 │   3000.0 │
│     2 │ 3000.00 │ 1000.00 │ 2000.00 │   3000.0 │
│     3 │ 3500.00 │       ∞ │       ∞ │   3000.0 │
└───────┴─────────┴─────────┴─────────┴──────────┘
```
"""
mutable struct LWModel{
    T1 <: AbstractArray{<:Any}, T2 <: AbstractArray{<:Any}, T3 <: AbstractArray{<:Any}} <:
               AbstractGeophyModel{T1, T2}
    m::T1
    h::T2
    ρ::T3
end

"""
create a model type for a given resistivity distribution
that can be used to calculate forward response for 1d MT

## Usage

```jldoctest
m = [4e3, 3e3, 3.5e3]
h = [1000.0, 1000.0]
d = [3e3, 3e3, 3e3]
vp = m .* 1.72
model = RWModel(m, h, d, vp)

# output

1D Rayleigh Wave Model
┌───────┬─────────┬─────────┬─────────┬──────────┬────────┐
│ Layer │      vₛ │       h │       z │        ρ │     vₚ │
│       │  [km/s] │     [m] │     [m] │ [kg⋅m⁻³] │ [km/s] │
├───────┼─────────┼─────────┼─────────┼──────────┼────────┤
│     1 │ 4000.00 │ 1000.00 │ 1000.00 │   3000.0 │ 6880.0 │
│     2 │ 3000.00 │ 1000.00 │ 2000.00 │   3000.0 │ 5160.0 │
│     3 │ 3500.00 │       ∞ │       ∞ │   3000.0 │ 6020.0 │
└───────┴─────────┴─────────┴─────────┴──────────┴────────┘
```
"""
mutable struct RWModel{T1 <: AbstractArray{<:Any}, T2 <: AbstractArray{<:Any},
    T3 <: AbstractArray{<:Any}, T4 <: AbstractArray{<:Any}} <: AbstractGeophyModel{T1, T2}
    m::T1
    h::T2
    ρ::T3
    vp::T4
end

RWModel(m, h, ρ, γ::T) where {T <: Number} = RWModel(m, h, ρ, m .* γ)

function SubsurfaceCore.default_params(::Type{T}) where {T <: Union{LWModel, RWModel}}
    (; mode=0, dc=0.01, dt=0.001, type=Val(:phase))
end
const default_params_surface_waves = default_params(RWModel) # same for LWModel
