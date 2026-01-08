```@raw html
---
# https://vitepress.dev/reference/default-theme-home-page
layout: home

hero:
  name: ProEM.jl
  text: Geophysical modeling and inversion in Julia
  tagline: Automatic Differentiation enabled probabilistic inference of subsurface
  image:
    src: logo.png
    alt: ProEM.jl
  actions:
    - theme: alt
      text: View on Github
      link: https://github.com/ayushinav/ProEM.jl
    
features:
  - icon: 🔢
    title: Modeling
    details: Estimate geophysical observables using rock physics
    link: intro/getting_started

  - icon: 📊
    title: Probabilstic inference
    details: Perform probabilistic inference of parameters
    link: /tutorials/stochastic_inverse

  - icon: ∂
    title: Differentiability
    details: Get derivatives using automatic differntiation
    link: /tutorials/ad
---
```

`ProEM.jl` is supposed to be a high performance code for doing forward and inverse modeling in geophysics using julia. We hope to write the package such that any other geophysical survey can also be used and we can tend towards a joint forward and inverse modeling library.

## Installation

You can install `ProEM.jl` on Julia by running:

```julia
using Pkg
Pkg.add("ProEM.jl")
```

## Available models

### Conductivity models

  - Magnetotellurics (MT)
  - DC resistivity

### Seismic models

  - Rayleigh waves
  - Love waves

TODO : Links to the models?
TODO : present probabilistic schemes here?
