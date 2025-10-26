function dnka!(C, wvno2, gam, gammk, rho, a0, cpcq, cpy, cpz, cqw, cqx, xy, xz, wy, wz)
    # constants
    gamm1 = gam - 1
    twgm1 = gam + gamm1
    gmgmk = gam * gammk
    gmgm1 = gam * gamm1
    gm1sq = gamm1 * gamm1

    a0pq = a0 -cpcq
    t = -2*wvno2

    C[1, 1] = cpcq - 2 * gmgm1 * a0pq - gmgmk * xz - wvno2 * gm1sq * wy
    C[1, 2] = (wvno2 * cpy - cqx) / rho
    C[1, 3] = -(twgm1 * a0pq + gammk * xz + wvno2 * gamm1 * wy) / rho
    C[1, 4] = (cpz - wvno2 * cqw) / rho
    C[1, 5] = -(2 * wvno2 * a0pq + xz + wvno2 * wvno2 * wy) / (rho*rho)

    C[2, 1] = (gmgmk * cpz - gm1sq * cqw) * rho
    C[2, 2] = cpcq
    C[2, 3] = gammk * cpz - gamm1 * cqw
    C[2, 4] = -wz
    C[2, 5] = C[1, 4]

    C[4, 1] = (gm1sq * cpy - gmgmk * cqx) * rho
    C[4, 2] = -xy
    C[4, 3] = gamm1 * cpy - gammk * cqx
    C[4, 4] = C[2, 2]
    C[4, 5] = C[1, 2]

    C[5, 1] = (
        -(2 * gmgmk * gm1sq * a0pq + gmgmk * gmgmk * xz + gm1sq * gm1sq * wy) * (rho*rho)
    )
    C[5, 2] = C[4, 1]
    C[5, 3] = (
        -(gammk * gamm1 * twgm1 * a0pq + gam * gammk * gammk * xz + gamm1 * gm1sq * wy)
        * rho
    )
    C[5, 4] = C[2, 1]
    C[5, 5] = C[1, 1]

    C[3, 1] = t * C[5, 3]
    C[3, 2] = t * C[4, 3]
    C[3, 3] = a0 + 2 * (cpcq - C[1, 1])
    C[3, 4] = t * C[2, 3]
    C[3, 5] = t * C[1, 3]

    return nothing;

end

function var(p, q, ra, rb, wvno, xka, xkb, dpth)
    pex = zero(p) # TODO
    if(wvno < xka)
        sinp = sin(p)
        w = sinp / ra
        x = -ra * sinp
        cosp = cos(p)
    elseif(wvno == xka)
        cosp = zero(p) + 1
        w = dpth
        x = zero(ra)
    elseif(wvno > xka)
        pex = p
        fac = exp(-2p) #ifelse(p < 16, exp(-2p), 0)
        cosp = (1 + fac) * 0.5f0
        sinp = (1 - fac) * 0.5f0
        w = sinp / ra
        x = ra * sinp
    end

    # Examine S-wave eigenfunctions
    # Checking whether c > vs, c = vs or c < vs
    sex = zero(q)
    if(wvno < xkb)
        sinq = sin(q)
        y = sinq / rb
        z = -rb * sinq
        cosq = cos(q)
    elseif(wvno == xkb)
        cosq = zero(q) + 1
        y = dpth
        z = zero(ra)
    elseif(wvno > xkb)
        sex = q
        fac = exp(-2q) #ifelse(q < 16, exp(-2q), 0)
        cosq = (1 + fac) * 0.5f0
        sinq = (1 - fac) * 0.5f0
        y = sinq / rb
        z = rb * sinq
    end

    # Form eigenfunction products for use with compound matrices
    exa = pex + sex
    a0 = exp(-exa) #ifelse(exa < 60, exp(-exa), zero(exa))
    cpcq = cosp * cosq
    cpy = cosp * y
    cpz = cosp * z
    cqw = cosq * w
    cqx = cosq * x
    xy = x * y
    xz = x * z
    wy = w * y
    wz = w * z

    qmp = sex - pex
    fac = exp(qmp) #ifelse(exa > -40, exp(qmp), zero(qmp))
    cosq *= fac
    y *= fac
    z *= fac

    return w, cosp, a0, cpcq, cpy, cpz, cqw, cqx, xy, xz, wy, wz
end

function dltar(c, ω, model::LWModel)
    vs = model.m
    h = model.h
    ρ = model.ρ
    
    xkb = ω / vs[end]

    rb = sqrt((c + xkb) * abs(c - xkb))
    e1 = ρ[end] * rb
    e2 = inv(vs[end]^2)

    for m in range(length(vs) - 1, 1, step = -1)
        
        G = ρ[m] * vs[m]^2
        xkb = ω / vs[m]

        rb = sqrt((c + xkb) * abs(c - xkb))
        q = h[m] * rb 

        if(c < xkb)
            sinq = sin(q)
            y = sinq / rb
            z = -rb * sinq
            cosq = cos(q)
        elseif(c == xkb)
            cosq = zero(q) + 1
            y = h[m]
            z = zero(q)
        else
            fac = ifelse(q < 16, exp(-2q), zero(q))
            cosq = (1 + fac) * 0.5f0
            sinq = (1 - fac) * 0.5f0
            y = sinq / rb
            z = rb * sinq
        end

        e10 = e1 * cosq + e2 * G * z
        e20 = e1 * y / G + e2 * cosq
        max_ = maximum([abs(e10), abs(e20)])
        dr = ifelse(max_ < 1f-40, zero(max_) + 1, max_)
        e1 = e10 / dr
        e2 = e20 / dr
    end

    return e1
end


function dltar(k, ω, model::RWModel)
    vp = model.vp
    vs = model.m
    ρ = model.ρ
    h = model.h

    e = zeros(1,5) # can be preallocated
    ee = zeros(1,5) # can be preallocated
    C = zeros(5,5) # can be preallocated

    (ω < 1f-4) && (ω = 1f-4 + zero(ω))

    xka = ω / vp[end]
    xkb = ω / vs[end]
    ra = sqrt((k + xka) * abs(k - xka))
    rb = sqrt((k + xkb) * abs(k - xkb))
    
    # t_ = vs[end] / ω

    # E matrix for the bottom half-space
    gammk = 2 * (vs[end] / ω)^2
    gam = gammk * k*k
    gamm1 = gam - 1
    ρ_end = ρ[end]

    e[1] = ρ_end^2 * (gamm1 * gamm1 - gam * gammk * ra * rb)
    e[2] = -ρ_end * ra
    e[3] = ρ_end * (gamm1 - gammk * ra * rb)
    e[4] = ρ_end * rb
    e[5] = k^2 - ra * rb

    # Matrix multiplication from bottom layer upward
    for m in range(length(vs) - 1, 1, step = -1)
        xka = ω / vp[m]
        xkb = ω / vs[m]
        # t = vs[m] / ω
        gammk = 2 * (vs[m] / ω)^2
        gam = gammk * k*k
        ra = sqrt((k + xka) * abs(k - xka))
        rb = sqrt((k + xkb) * abs(k - xkb))

        dpth = h[m] # should change later on
        p = ra * dpth
        q = rb * dpth

        # Evaluate cosP, cosQ...
        _, _, a0, cpcq, cpy, cpz, cqw, cqx, xy, xz, wy, wz = var(
            p, q, ra, rb, k, xka, xkb, dpth
        )

        # Evaluate Dunkin's matrix
        dnka!(C,
            k*k, gam, gammk, ρ[m], a0, cpcq, cpy, cpz, cqw, cqx, xy, xz, wy, wz
        )

        mul!(ee, e, C);
        norm_fac = maximum(abs.(ee))
        e .= ee./norm_fac


    end

    # xka = ω / vp[1]
    # ra = sqrt((c + xka) * abs(c - xka))
    # dpth = h[1]
    # p = ra * dpth

    # rb = sqrt((c + xkb) * abs(c - xkb))
    # q = rb * dpth
    # w, cosp, _, _, _, _, _, _, _, _, _, _ = var(
    #     p, q, ra, 1f-5, c, xka, xkb, dpth
    # )
    
    return e[1]

end

function get_c(t, m, mode, dc)

    c = zero(t)
    c = zero(t)

    b_range = extrema(m.m)
    c_high = 10 # should not go this far really
    c_low = b_range[1] .* 0.8
    
    resp_ = zero(t)
    
    err_fn(c, omega) = dltar(omega / c, omega, m)

    for i in eachindex(t) # this can be parallelized
        
        ω = 2π/t[i]
        f(c, p) = dltar(ω/c, ω, m) #err_fn(c, 2π/t[i])
        
        c_low = b_range[1] .* 0.8
        c_high_each = c_low

        for im in 1:mode+1

            f_low = f(c_low, [])
            while c_high_each <= c_high
                f_high_each = f(c_high_each, [])
                if sign(f_high_each) * sign(f_low) < 0
                    break
                else
                    c_high_each += dc
                end

                if c_high_each> c_high
                    @warn "search space exceeded! t = $(t[i])"
                    break;
                end
            end

            prob_init = IntervalNonlinearProblem(f, (c_high_each - 2dc, c_high_each), [])
            sol = solve(prob_init)
            c = sol.u
            
            c_low  = c + dc
            resp_[i] = c

        end
        
    end
    return resp_

end