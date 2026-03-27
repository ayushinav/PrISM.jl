# Stochastic inversion

Stochastic inversion revolves around the Bayesian formulation :

```math
p(m|d) \propto p(d|m) p(m)
```

where $d$ is the data to be inverted and $m$ is the model. $p(m)$ defines the prior information about the model (our guess about what we know about the model prior to performing inversion). We update this to obtain the posterior distribution $p(m|d)$ using the observed dataset. $p(d|m)$ is called the likelihood term and includes the physics of the system as well as the metric of misfit.

We do not generally get an explicit analytical form of the posterior distribution, but obtain the samples from it. The most common way to obtain these samples is via Markov chain Monte Carlo (MCMC). Without going into details of MCMC, we complete here by mentioning that we usually need a large number of samples to obtain a good estimate of the posterior distribution. Naive MCMC schemes such as Random Walk Metropolis can take a long time to converge, and using gradient information can lead to a more efficient traversal of the probability space.

To sum it up, performing stochastic inversion involves

  - specifying the *a priori* distribution
  - specifying the likelihood (physics + misfit)
  - specifying the MCMC scheme (sampler + no of samples)

## Constructing distributions

### Model distribution (*a priori* information)

Before beginning to talk about how to construct the *a priori* distribution, it is important to understand that the model here refers to the discretization space as well as the values of physical properties. For e.g., for electrical methods, the model space will consist of the electrical resistivity as well as the grid sizes discretizing the spatial domain. You can vary both or keep the grid size fixed. Grid sizes are usually kept fixed. Varying them can be useful in 1D, but it can also have its implicit effect on the results

For most applications, however, we fix the node points. We follow a 1D MT example to show the framework. For now, we begin by choosing a very broad prior for all the `n` layers. A model distribution can be constructed using `MTModelDistribution(...)`. This has the same structure as `MTModel`, that is, the first parameter denotes the prior for electrical conductivities, while the second is for the layer thicknesses `h`. If you do not want to infer on `h`, just pass it as a simple vector, as is also demonstrated in the following case:

```julia
n = 50 # number of layers
h = fill(20.0, n - 1) # making the discretization

# prior space with fixed model discretization
modelD = MTModelDistribution(Product([Uniform(-1, 5) for i in 1:n]), vec(h));
```

Some of the things to understand here are :

  - `Product(...)` is a `Distributions.jl` function that joins a series of multivariate/univariate samples, here `Uniform(-1,5)`. Thus, `Product([Uniform(-1, 5) for in 1:n])` makes a sampler that will output an `n`-length vector, with each element sampled from `Uniform(-1,5)`. Here, we use a uniform prior, but one can choose any prior. The only thing to keep in mind is that the prior can be sampled, that is, `rand(modelD.m)` returns a vector of appropriate length. Another eg., would be `MultivariateNormal(μ, Σ)`.
  - To also sample `h`, pass another multivariate distribution that outputs `n-1` length vector. Again using a uniform, independent prior, we'll have:

```julia

# prior space with variable model discretization
modelD = DCModelDistribution(
    Product([Uniform(-1, 5) for i in 1:n]), Product([Uniform(15, 25) for i in 1:(n - 1)]));
```

The above defines the *a priori* distribution for DC resistivity modeling. Corresponding functions exist for other geophysical models.

### Response distribution (*likelihood*)

A likelihood is determined by an observed response `r_obs` we want to fit, and the errors associated with it `err_resp`. One of the popular ways in which likelihood can be formed is using the gaussian distribution, centered around `r_obs` with variance given by `err_resp`. This is the probabilistic equivalence of mean squared misfit.

Remember that to make a likelihood, we need the physics of the system and the error metric. With the data at avail, we now need to specify the misfit function. Here, we just need to pass a function that can take in a response parameter and the associated error and produce a distribution. We provide a function `norm_dist` that takes in a vector for the response and another vector/matrix for the covariance matrix of the errors. The response distribution is then constructed by:

```julia
respD = MTResponseDistribution(normal_dist, normal_dist)
```

!!! note
    
    Both `r_obs` and `err_resp` have the same type, eg. `MTResponse`.

In the above, `MTResponseDistribution` specifies the physics of the system to completely specify the likelihood. Similar `ResponseDistribution`s exist for other type of geophysical models.

## Inference

We now have all the ingredients to perform inversion, except the sampler. This [page](https://turinglang.org/docs/usage/sampler-visualisation/) provides a brief review of samplers `Turing.jl` provides. The number of posterior points to be sampled `n_samples`, the algorithm `mcmc_alg` and the distributions are brought together by `mcmc_cache`.

```julia
mcmc_alg = NUTS();
n_samples = 10_000;
mcache = mcmc_cache(modelD, respD, n_samples, mcmc_alg);
```

The posterior samples are then sampled by simply calling:

```julia
mcmc_chain = stochastic_inverse(r_obs, err_resp, ω, mcache)
```
