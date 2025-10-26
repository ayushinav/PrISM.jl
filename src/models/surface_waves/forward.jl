
const default_surf_tf_fns = (;c = lin_tf)

function surf96(t, m, mode, dc, dt, ::Val{:phase})
    return get_c(t, m, mode, dc)
end

function surf96(t, m, mode, dc, dt, ::Val{:group})
    @. t = t * inv(1 + dt)
    c_plus = get_c(t, m, mode, dc)

    @. t = t * inv(1 - dt)
    c1_minus = get_c(t, m, mode, dc)

    c = zeros(typeof(t[1] * m.m[1]), length(t))
    for i in eachindex(t)
        
    end

    c = @. 2dt * inv(t) * ((1+dt)* inv(t * c_plus) - (1-dt)* inv(t * c_minus))
    return c

end

function SubsurfaceCore.forward(m::Tm, t::T3, response_trans_utils::T=default_surf_tf_fns, params=default_params_surface_waves) where {Tm <: Union{LWModel, RWModel}, T, T3}
    c = surf96(t, m, params.mode, params.dc, params.dt, Val(params.type))
    SurfaceWaveResponse(response_trans_utils.c.tf.(c))
end
