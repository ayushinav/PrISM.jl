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

function SubsurfaceCore.forward(m::Tm, t::T3, response_trans_utils::T=default_surf_tf_fns,
        params=default_params_surface_waves) where {Tm <: LWModel, T, T3}
    c = zeros(eltype(m.m), length(t))
    surf96!(c, t, m, params.mode, params.dc, params.dt, Val(params.type))
    f1 = response_trans_utils.c.tf
    SurfaceWaveResponse{typeof(c)}(f1.(c))
end

function SubsurfaceCore.forward(m::Tm, t::T3, response_trans_utils::T=default_surf_tf_fns,
        params=default_params_surface_waves) where {Tm <: RWModel, T, T3}
    c = zeros(eltype(m.m), length(t))
    surf96!(c, t, m, params.mode, params.dc, params.dt, Val(params.type))
    f1 = response_trans_utils.c.tf
    SurfaceWaveResponse{typeof(c)}(f1.(c))
end

function forward!(resp::Tr,
        m::Tm,
        t::T3,
        response_trans_utils::T=default_surf_tf_fns,
        params=default_params_surface_waves) where {
        Tm <: RWModel, T, T3, Tr <: SurfaceWaveResponse}
    surf96!(resp.c, t, m, params.mode, params.dc, params.dt, Val(params.type))
    f1 = response_trans_utils.c.tf
    broadcast!(f1, resp.c, resp.c)
    return nothing
end

function forward!(resp::Tr,
        m::Tm,
        t::T3,
        response_trans_utils::T=default_surf_tf_fns,
        params=default_params_surface_waves) where {
        Tm <: LWModel, T, T3, Tr <: SurfaceWaveResponse}
    surf96!(resp.c, t, m, params.mode, params.dc, params.dt, Val(params.type))
    f1 = response_trans_utils.c.tf
    broadcast!(f1, resp.c, resp.c)
    return nothing
end
