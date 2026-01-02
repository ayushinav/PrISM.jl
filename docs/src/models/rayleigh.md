
# Rayleigh wave (MT)

```@setup rw_demo
using ProEM, CairoMakie, InteractiveUtils
```

## Model

We assume the following subsurface resistivity distribution with 4 layers:
6.8  5.95  7.14  7.65
| Layer # | thickness (m) | $V_s \; (km s^{-1}$) | $V_p \; (km s^{-1}$) | $\rho \; (g cm^{-3}$) | 
| :-: | :-: | :-: | :-: | :-: |
| 1 | 8000 | 4 | 6.8 | 2
| 2 | 4000 | 3.5 | 5.95 | 2
| 3 | 8000 | 4.2 | 7.14 | 2
| 4 | $\infty$ | 4.5 | 7.65 | 2


The `RWModel` is defined as :

```@example rw_demo
vs = [4, 3.5, 4.2, 4.5]
vp = vs .* 1.7
ρ = [2., 2., 2., 2.]
h = [8e3, 4e3, 8e3]
m = RWModel(vs, h, ρ, vp)
```

```@raw html
<details closed><summary>Code for this figure</summary>
```

```@example rw_demo
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

To obtain the responses, that is apparent resistivity `ρₐ` and phase `ϕ`, we first need to define the frequencies :

```@example rw_demo
freq = exp10.(-1:0.2:2)
T = inv.(freq)
```

and then call the `forward` dispatch and get `MTResponse`:

```@example rw_demo
resp = forward(m, T)
nothing # hide
```

```@raw html
<details closed><summary>Code for this figure</summary>
```

```@example rw_demo
fig = Figure() #size = (800, 600))
ax = Axis(fig[1, 1]; backgroundcolor=(:magenta, 0.05), aspect=AxisAspect(1))

plot_response!([ax], T, resp; plt_type=:scatter)
nothing # hide
```

```@raw html
</details>
```

```@example rw_demo
fig # hide
```

### In-place operations

Mutating forward calls are also supported. This leads to no extra allocations while calculating the forward responses. This also speeds up the performance, though the calculations in the well-optimized forward calls are way more expensive than allocations to get a significant boost here. We now call the in-place variant `forward!` and provide an `MTResponse` variable to be overwritten :

```@example rw_demo
forward!(resp, m, T)
nothing # hide
```

We can compare the allocations made in the two calls:

```@example rw_demo
@allocated forward!(resp, m, T)
```

```@example rw_demo
@allocated forward(m, T)
```

## Benchmarks

While the exact runtimes will vary across processors, the runtimes increase linearly with number of layers.

```@raw html
<details closed><summary>Benchmark code</summary>
```

```@example rw_demo
n_layers = Int.(exp2.(2:1:6))
freq = exp10.(-1:0.2:2)
T = 2π .* freq

btime_iip = Float64.(zero(n_layers))
btime_oop = Float64.(zero(n_layers))
for i in eachindex(n_layers)
    vs = 3 .+ 0.5 .* randn(n_layers[i])
    vp = vs .* 1.72
    ρ = 2 .* ones(n_layers[i])
    h = 100 .* rand(n_layers[i]-1)
    m = RWModel(vs, h, vp, ρ)
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

```@example rw_demo
fig # hide
```

### Reproducibility

```@example rw_demo
println("The above benchmarks were obtained on $(Sys.cpu_info()[1].model)") # hide
# println(versionfo()) # hide
```

Hello Mehdia,
Happy birthday and a very happy new year! I have always admired you for the person you are. Honestly, I consider myself lucky to have found a friend in you. I think this year was a great milestone for you and a testament to your own commitments. You passed your comprehensive exam, not that there was any doubt, but it is still a hard objective you have achieved since coming here. I have already expressed my awe at the transformation because of working out, but you have achieved even greater transformation in other aspects. That picture you showed us in front of your school, with 10 years in between, captures beautifully how long you have come. I have always appreciated our discussions on relationships and religion, and I hope we continue to do so. In a world getting more and more radicalised, you are a gem to talk critically about the values you have held close growing up. I hope the coming year brings you joy, and fills your heart with love and happiness. 