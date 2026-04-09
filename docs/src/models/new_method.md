# Adding a new method

If you want to add a new method, please open an issue/ pull request.
We demonstrate adding a new method called `MyMethod` which depends on the physical property `myProperty`, and generates the response `myResponse`. For the sake of brevity, we assume that `MyProperty` is the only parameter we are interested in. Certain models such as surface waves depend on shear-wave velocity along with p-wave velocity and density. Similarly, magnetotelluric method has two responses, apparent resistivity and phase.

!!!tip
    The easiest way to understand this is by going through the following steps along with the model definitions of `MTModel` and `MTResponse`.

1. Create a new file `models/myproperty/types.jl` and define 

```julia
mutable struct MyProperty{T1, T2} <: AbstractGeophyModel
    m::T1 # m is the parameter corresponding to the physical property
    h::T2 # h is the discretization, e.g. layer thicknesses
end
```

for the geophysical model, and similarly 

```julia
mutable struct MyResponse{T1, T2} <: AbstractGeophyResponse
    r::T1
end
```

for the corresponding response. The variables should be strictly named `m` and `h` for the geophysical model, whereas, they can be anything for the response, here `r`.

2. Define defaults for modeling. This could be anything such as discretization in solving a PDE related to the system, step size in a grid search algorithm or a threshold to achieve to assume convergence of a system

```julia
SubsurfaceCore.default_params(::Type{T}) where {T <: MyModel} = (; stepsize = 2., threshold = 1e-8)
const default_params_mymodel = default_params(MyModel)
```

3. Create a new file `models/forward.jl`
These defaults are then used to calculate the forward response of the system at say, certain frequencies. We define two versions :
**Non-mutating functions** : they return the response, e.g.

```julia
function SubsurfaceCore.forward(m, freqs, params = default_params)
    # some computation that calculates the response vector, say resp
    return MyResponse(resp)
end
```

**Mutating functions ** : they take the response object and mutate it, e.g.

```julia
function SubsurfaceCore.forward!(my_resp, m, freqs, params = default_params)
    # some computation that modifies the my_resp structure
    return nothing
end
```

4. Define the default transformation functions. This is particularly useful for Bayesian inference when one of your data is on, say log-scale.

5. Define a `mutable struct` ` MyModelDistribution` under `probabilistic/models` with the same form as for step 2.

