# Occam1D

```@setup occam_demo
using ProEM, LinearAlgebra, CairoMakie
```
## Brief introduction

We provide inversion algorithm popularized as [Occam1D](https://marineemlab.ucsd.edu/steve/bio/Occam1D.pdf).

Occam1D tries to find the model with smoothest variations that fits the data. It does so by trying to minimize the first derivative of the model vector along with the misfit function. Unsurprisingly, this is done by modifying the regularizer term as

```math
\begin{align*}
& m^T C_m m \\
= & m^T (L^T L) m \\
= & (Lm)^T (Lm) \\
=& || Lm ||_2^2
\end{align*}
```

Occam1D also iterates over $\mu$ at every step to find the largest $\mu$ that best fits the data, steering it towards smooth models. Therefore, in every iteration, it solves for a linearized solution as well as the regularization coefficient.

## Demo

We demo using Occam on a couple of models below:

### Magnetotellurics (MT)

Occam1D was initially designed and compiled for 1D MT, and this is the first example we run here.

Let's define a synthetic model:

```@example occam_demo
ρ = log10.([100.0, 10.0, 400.0, 1000.0])
h = [100.0, 1000.0, 10000.0]
m = MTModel(ρ, h)
```

and then get data, and assume 5% error floors:

```@example occam_demo
T = 10 .^ (range(-1, 5; length=19))
ω = 2π ./ T
resp = forward(m, ω)

err_resp = MTResponse(
    0.1 .* resp.ρₐ,
    180 / π .* asin(0.1) .+ zero(ω)
)
```

We need to define a data covariance matrix $C_d$ and an initial model. Let's assume gaussian noise in this case, and choose an initial model a half-space of 100 $\Omega m$.

```@example occam_demo
C_d = diagm(inv.([err_resp.ρₐ..., err_resp.ϕ...])) .^ 2; # weight matrix

# initial model
h_test = 10 .^ range(0.0, 5.0; length=50)
ρ_test = 2 .* ones(length(h_test) + 1)

m_occam = MTModel(ρ_test, h_test);
```

The final result will be stored in the same variable `m_occam`. All we need to do now is specify using Occam and then calling `inverse!`.

!!! note
    Using Occam, you also have the option to perform a smoothing step. Once the model has achieved the threshold misfit, it is smoothened until it fits the data just about the threshold misfit.

```@example occam_demo
alg_cache = Occam(; μgrid=[1e-2, 1e6])
inverse!(m_occam, resp, ω, alg_cache; W=C_d, max_iters=50, verbose=true, smoothing_step = true)
```

We can then plot the data and the models to see the fits.


```@example occam_demo
fig = Figure()
ax_m = Axis(fig[:,1])

plot_model!(ax_m, m, label = "true", color = :steelblue3)
plot_model!(ax_m, m_occam, label = "occam", color = :tomato)


ax1 = Axis(fig[1,2])
ax2 = Axis(fig[2,2])

plot_response!([ax1, ax2], T, resp; plt_type = :scatter, color = :steelbue3, label = "true")
plot_response!([ax1, ax2], T, resp; errs = err_resp, plt_type = :errors, whiskerwidth=10, color = :steelbue3)

resp_occam = forward(m_occam, ω)
plot_response!([ax1, ax2], T, resp; color = :tomato, label = "occam")
```