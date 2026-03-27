# Vizualization

## Model

You can plot models using the dispatchs

```@docs; cannonical = false
plot_model
plot_model!
```

!!! note
    
    Do note that the above functions only plot the parameter denoted by `m`, even though there might be other parameters, e.g. for Rayleigh wave models, we invert of shear wave velocity even though the model also requires p-wave velocities and densities.

## Response

Similarly, to plot response for different geophysical surveys, we have

```@docs; cannonical = false
plot_response
plot_response!
```
