# DC (Direct Current) Resistivity

```@setup dc_demo
using PrISM, CairoMakie, InteractiveUtils
```

## Model

We assume the following subsurface resistivity distribution with 4 layers:

| Layer # | thickness (m) | $\rho \; (\Omega m$) |
|:-------:|:-------------:|:--------------------:|
| 1       | 100           | 500                  |
| 2       | 200           | 1000                 |
| 3       | 20            | 10                   |
| 4       | $\infty$      | 100                  |

The `DCModel` is defined as :

```@example dc_demo
ρ = log10.([500.0, 1000.0, 10.0, 100.0])
h = [100.0, 200.0, 20.0]
m = DCModel(ρ, h)
```

```@raw html
<details closed><summary>Code for this figure</summary>
```

```@example dc_demo
fig = Figure()
ax = Axis(fig[1, 1]; xlabel="log ρ (locsm)", ylabel="depth (m)", backgroundcolor=(
    :magenta, 0.05))

plot_model!(ax, m)
nothing # hide
```

```@raw html
</details>
```

```@example dc_demo
fig # hide
```

!!! warn
    
    Always use `Float64` or `Float32` types while defining the vectors for resistivity and thickness. Using `Int` will throw an `InexactError`, e.g. : `InexactError: Int64(4193.453970907305)`

## Response

To obtain the responses, that is apparent resistivity `ρₐ`, we first need to define the current source and electrodes (receivers) locations. We provide two most common source receiver configurations

```@docs; canonical = false
get_schlumberger_array
get_wenner_array
```

Here, we demonstrate using the Wenner array for a vertical sounding survey:

```@example dc_demo
locs = get_wenner_array(100:50:2000)
```

and then call the `forward` dispatch as usual to get `MTResponse`:

```@example dc_demo
resp = forward(m, locs)
```

```@raw html
<details closed><summary>Code for this figure</summary>
```

```@example dc_demo
fig = Figure() #size = (800, 600))
ax1 = Axis(fig[1, 1]; backgroundcolor=(:magenta, 0.05), aspect=AxisAspect(1))

ab_2 = abs.(locs.srcs[:, 2] .- locs.srcs[:, 1]) ./ 2
plot_response!([ax1], ab_2, resp; plt_type=:scatter)
nothing # hide
```

```@raw html
</details>
```

```@example dc_demo
fig # hide
```

### In-place operations

Mutating forward calls are also supported. This leads to no extra allocations while calculating the forward responses. This also speeds up the performance, though the calculations in the well-optimized forward calls are way more expensive than allocations to get a significant boost here. We now call the in-place variant `forward!` and provide an `MTResponse` variable to be overwritten :

```@example dc_demo
forward!(resp, m, locs)
nothing # hide
```

We can compare the allocations made in the two calls:

```@example dc_demo
@allocated forward!(resp, m, locs)
```

```@example dc_demo
@allocated forward(m, locs)
```

## Benchmarks

While the exact runtimes will vary across processors, the runtimes increase linearly with number of layers.

```@raw html
<details closed><summary>Benchmark code</summary>
```

```@example dc_demo
n_layers = Int.(exp2.(2:1:6))
locs = get_wenner_array(100:50:2000)

btime_iip = Float64.(zero(n_layers))
btime_oop = Float64.(zero(n_layers))
for i in eachindex(n_layers)
    ρ = 2 .* randn(n_layers[i])
    h = 100 .* rand(n_layers[i]-1)
    m = DCModel(ρ, h)
    resp = forward(m, locs)
    time_oop = @timed begin
        for i in 1:1000
            forward(m, locs)
        end
    end
    btime_oop[i] = time_oop.time/1e3
    forward!(resp, m, locs)
    time_iip = @timed begin
        for i in 1:1000
            forward!(resp, m, locs)
        end
    end
    btime_iip[i] = time_iip.time/1e3
end

fig = Figure()
ax = Axis(fig[1, 1]; xscale=log2, backgroundcolor=(:magenta, 0.05),
    xlabel="no. of layers", ylabel="time (ms)")
lines!(n_layers, btime_oop .* 1e3; label="out-of place forward calls", color=:tomato)
lines!(n_layers, btime_iip .* 1e3; label="in place forward calls",
    linestyle=:dash, color=:steelblue3)
Legend(fig[1, 2], ax)
nothing # hide
```

```@raw html
</details>
```

```@example dc_demo
fig # hide
```

### Reproducibility

```@example dc_demo
println("The above benchmarks were obtained on $(Sys.cpu_info()[1].model)") # hide
# println(versionfo()) # hide
```
