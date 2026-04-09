```@raw html
---
# https://vitepress.dev/reference/default-theme-home-page
layout: home

hero:
  name: PrISM.jl
  text: PRobabilistic Inference of Subsurface Models in Julia
  tagline: Automatic Differentiation enabled deterministic and stochastic imaging of subsurface using geophysical data
  image:
    src: logo.png
    alt: PrISM.jl
  actions:
    - theme: alt
      text: View on Github
      link: https://github.com/ayushinav/PrISM.jl
    
features:
  - icon: 🔢
    title: Modeling
    details: Estimate subsurface models using deterministic inversion using wide suite of solvers
    link: deterministic_inverse/index

  - icon: 📊
    title: Probabilstic inference
    details: Perform Uncertainty Quantification of the subsurface models
    link: stochastic_inverse/index

  - icon: ∂
    title: Differentiability
    details: Get derivatives using automatic differntiation
    link: /tutorials/ad
---
```

`PrISM.jl` is a high performance code for doing forward and inverse modeling in geophysics using julia. While the current capabilities sit at 1D methods, we hope to extend it for 2D and 3D problems as well.

## Installation

You can install `PrISM.jl` on Julia by running:

```julia
using Pkg
Pkg.add("PrISM.jl")
```

## Available models

### Conductivity models

  - Magnetotellurics (MT)
  - DC resistivity

### Seismic models

  - Rayleigh waves
  - Love waves

