const default_surf_tf_fns = (; c=no_tf)

function surf96!(c, t, m, mode, dc, dt, ::Val{:phase})
    rmul!(m.h, inv(eltype(m.h)(1000)))
    get_c!(c, t, m, mode, dc)
    rmul!(m.h, eltype(m.h)(1000))
    return nothing
end

function surf96!(c, t, m, mode, dc, dt, ::Val{:group}) # TODO
    rmul!(m.h, inv(eltype(m.h)(1000)))
    @. t = t * inv(1 + dt)
    c_plus = zero(c)
    get_c!(c_plus, t, m, mode, dc)

    @. t = t * inv(1 - dt)
    get_c!(c, t, m, mode, dc) # c_minus 

    @. c = 2dt * inv(t) * ((1 + dt) * inv(t * c_plus) - (1 - dt) * inv(t * c))
    rmul!(m.h, eltype(m.h)(1000))
    return nothing
end

"""
   forward(m::LWModel, t)

returns `LWResponse` for the given `LWModel m` at the periods defined by `t`

## Arguments

 - `m` : LWModel for forward response
 - `t` : periods for which forward response is to be calculated

## Optional arguments
 - `params` : `NamedTuple` containing configuration for forward calculation, contains 
    - `mode` : dispersion curve mode, defaults to `0` implying fundamental node
    - `dc` : step size used to search for grid space solution of propagator matrix, defaults to `0.001`
    - `dt` : Δt to be used for calculation of group velocity from phase velocity
    - `type` : `Val` type variable that specifies whether to calculate group velocity or phase. Available options are :
        - `Val(:phase)` : phase velocity (default)
        - `Val(:group)` : group velocity
"""
function SubsurfaceCore.forward(m::Tm, t::T3, params=default_params_surface_waves) where {Tm <: LWModel, T3}
    c = zeros(eltype(m.m), length(t))
    surf96!(c, t, m, params.mode, params.dc, params.dt, params.type)
    SurfaceWaveResponse(c)
end

"""
   forward(m::RWModel, t)

returns `RWResponse` for the given `RWModel m` at the periods defined by `t`

## Arguments

 - `m` : RWModel for forward response
 - `t` : periods for which forward response is to be calculated

## Optional arguments
 - `params` : `NamedTuple` containing configuration for forward calculation, contains 
    - `mode` : dispersion curve mode, defaults to `0` implying fundamental node
    - `dc` : step size used to search for grid space solution of propagator matrix, defaults to `0.001`
    - `dt` : Δt to be used for calculation of group velocity from phase velocity
    - `type` : `Val` type variable that specifies whether to calculate group velocity or phase. Available options are :
        - `Val(:phase)` : phase velocity (default)
        - `Val(:group)` : group velocity
"""
function SubsurfaceCore.forward(m::Tm, t::T3, params=default_params_surface_waves) where {Tm <: RWModel, T3}
    c = zeros(eltype(m.m), length(t))
    surf96!(c, t, m, params.mode, params.dc, params.dt, params.type)
    SurfaceWaveResponse(c)
end

"""
   forward!(resp::LWResponse, m::LWModel, t)

overwrites `LWResponse` for the given `LWModel m` at the periods defined by `t`

## Arguments

 - `resp` : `LWResponse` to be overwritten
 - `m` : LWModel for forward response
 - `t` : periods for which forward response is to be calculated

## Optional arguments
 - `params` : `NamedTuple` containing configuration for forward calculation, contains 
    - `mode` : dispersion curve mode, defaults to `0` implying fundamental node
    - `dc` : step size used to search for grid space solution of propagator matrix, defaults to `0.001`
    - `dt` : Δt to be used for calculation of group velocity from phase velocity
    - `type` : `Val` type variable that specifies whether to calculate group velocity or phase. Available options are :
        - `Val(:phase)` : phase velocity (default)
        - `Val(:group)` : group velocity
"""
function forward!(resp::Tr,
        m::Tm,
        t::T3,
        params=default_params_surface_waves) where {
        Tm <: LWModel, T3, Tr <: SurfaceWaveResponse}
    surf96!(resp.c, t, m, params.mode, params.dc, params.dt, params.type)
    return nothing
end

"""
   forward!(resp::RWResponse, m::RWModel, t)

overwrites `RWResponse` for the given `RWModel m` at the periods defined by `t`

## Arguments

 - `resp` : `RWResponse` to be overwritten
 - `m` : LWModel for forward response
 - `t` : periods for which forward response is to be calculated

## Optional arguments
 - `params` : `NamedTuple` containing configuration for forward calculation, contains 
    - `mode` : dispersion curve mode, defaults to `0` implying fundamental node
    - `dc` : step size used to search for grid space solution of propagator matrix, defaults to `0.001`
    - `dt` : Δt to be used for calculation of group velocity from phase velocity
    - `type` : `Val` type variable that specifies whether to calculate group velocity or phase. Available options are :
        - `Val(:phase)` : phase velocity (default)
        - `Val(:group)` : group velocity
"""
function forward!(resp::Tr,
        m::Tm,
        t::T3,
        params=default_params_surface_waves) where {
        Tm <: RWModel, T3, Tr <: SurfaceWaveResponse}
    surf96!(resp.c, t, m, params.mode, params.dc, params.dt, params.type)
    return nothing
end
