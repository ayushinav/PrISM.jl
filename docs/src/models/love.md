

# Love wave

```@setup lw_demo
using ProEM, CairoMakie, InteractiveUtils
```

## Model

We assume the following subsurface resistivity distribution with 4 layers:

| Layer # | thickness (m) | $V_s \; (km s^{-1}$) | $\rho \; (g cm^{-3}$) | 
| :-: | :-: | :-: | :-: |
| 1 | 8000 | 4 | 2 |
| 2 | 4000 | 3.5 | 2 |
| 3 | 8000 | 4.2 | 2 |
| 4 | $\infty$ | 7.65 | 2 |


The `LWModel` is defined as :

```@example lw_demo
vs = [3.2, 4.39731, 4.40192, 4.40653, 4.41113, 4.41574]
ρ = [2.6, 3.38014, 3.37797, 3.37579, 3.37362, 3.37145]
h = [20.0, 20.0, 20.0, 20.0, 20.0] .* 1e3
m = LWModel(vs, h, ρ)
```

```@raw html
<details closed><summary>Code for this figure</summary>
```

```@example lw_demo
fig = Figure()
ax = Axis(fig[1, 1]; xlabel="Vs (km/s)", ylabel="depth (m)", backgroundcolor=(
    :magenta, 0.05))

plot_model!(ax, m)
nothing # hide
```

```@raw html
</details>
```

```
fig # hide
```

!!! warn
    
    Always use `Float64` or `Float32` types while defining the vectors for resistivity and thickness. Using `Int` will throw an `InexactError`, e.g. : `InexactError: Int64(4193.453970907305)`

## Response

To obtain the response, that is the dispersion curve `c`. We first need to define the frequencies :

```@example lw_demo
freq = exp10.(-2:0.1:1)
T = inv.(freq)
```

and then call the `forward` dispatch and get `SurfaceWaveResponse`:

```@example lw_demo
resp = forward(m, T)
nothing # hide
```

```@raw html
<details closed><summary>Code for this figure</summary>
```

```@example lw_demo
fig = Figure() #size = (800, 600))
ax = Axis(fig[1, 1]; backgroundcolor=(:magenta, 0.05), aspect=AxisAspect(1))

plot_response!([ax], T, resp; plt_type=:scatter)
nothing # hide
```

```@raw html
</details>
```

```@example lw_demo
fig # hide
```

!!! note
    By default, the forward response calculates the phase velocity at fundamental mode. This information is passed through params.
    ```@example rw_demo
    default_params(RWModel)
    ```

    Various symbols are defined as : 
    * `mode` : fundamental mode corresponds to 0, higher modes correspond to 1,2,...
    * `dc` : step size used to obtain solution of the propagater matrix
    * `dt` : step size used to obtain the derivative for group velocity, unused for phase velocity
    * `Val(:phase)` : to calculate phase velocity, pass `Val(:group)` to calculate group velocity

### In-place operations

Mutating forward calls are also supported. This leads to no extra allocations while calculating the forward responses. This also speeds up the performance, though the calculations in the well-optimized forward calls are way more expensive than allocations to get a significant boost here. We now call the in-place variant `forward!` and provide a `SurfaceWaveResponse` variable to be overwritten :

```@example lw_demo
forward!(resp, m, T)
nothing # hide
```

We can compare the allocations made in the two calls:

```@example lw_demo
@allocated forward!(resp, m, T)
```

```@example lw_demo
@allocated forward(m, T)
```

## Benchmarks

While the exact runtimes will vary across processors, the runtimes increase linearly with number of layers.

```@raw html
<details closed><summary>Benchmark code</summary>
```

```@example lw_demo
n_layers = Int.(exp2.(2:1:6))
freq = exp10.(-1:0.2:2)
T = 2π .* freq

btime_iip = Float64.(zero(n_layers))
btime_oop = Float64.(zero(n_layers))
for i in eachindex(n_layers)
    vs = 3 .+ 0.005 .* randn(n_layers[i])
    ρ = 2 .* ones(n_layers[i])
    h = 100 .* rand(n_layers[i]-1)
    m = LWModel(vs, h, ρ)
    resp = forward(m, T)
    time_oop = @timed begin 
        for i in 1:1000 forward(m, T) end
    end
    btime_oop[i] = time_oop.time/1e3
    forward!(resp, m, T)
    time_iip = @timed begin
        for i in 1:1000 forward!(resp, m, T) end
    end
    btime_iip[i] = time_iip.time/1e3
end

fig = Figure()
ax = Axis(fig[1,1], xscale = log2, backgroundcolor = (:magenta, 0.05), xlabel = "no. of layers", ylabel = "time (ms)")
lines!(n_layers, btime_oop .* 1e3, label = "out-of place forward calls", color = :tomato)
lines!(n_layers, btime_iip .* 1e3, label = "in place forward calls", linestyle = :dash, color = :steelblue3)
Legend(fig[1,2], ax)
nothing # hide
```

```@raw html
</details>
```

```@example lw_demo
fig # hide
```

### Reproducibility

```@example lw_demo
println("The above benchmarks were obtained on $(Sys.cpu_info()[1].model)") # hide
# println(versionfo()) # hide
```
