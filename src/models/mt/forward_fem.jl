function get_Z(m::MTModel{T1, T2}, ω) where {T1 <: AbstractVector, T2 <: GridapType}
    order = 1
    reffe = ReferenceFE(lagrangian, Float64, order)

    V = TestFESpace(
        m.h, reffe; dirichlet_tags=["tag1", "tag2"], vector_type=Vector{ComplexF64})

    U = TrialFESpace(V, [1 + 0im, 0 + 0im])

    Ω = Triangulation(m.h)
    dΩ = Measure(Ω, 2)

    f_val = 0.0 # RHS

    n_cells = num_cells(m.h)
    ρ_values_per_cell = zeros(Float64, n_cells)

    ρ_mat = reshape(view(ρ_values_per_cell, :), n_cells)
    ρ_mat .= m.m
    ρ_discrete = CellField(ρ_values_per_cell, Ω)

    # Setup
    a(u, v) = ∫(∇(u) ⋅ ∇(v) + im * ω * μ * inv(ρ_discrete) * u ⋅ v)dΩ
    b(v) = ∫(f_val * v)dΩ

    op = AffineFEOperator(a, b, U, V)
    E = Gridap.solve(op)

    E_interp = interpolate(E, V)

    E1 = E_interp(Fields.Point(-1.0))
    E2 = E_interp(Fields.Point(0.0))
    E3 = E_interp(Fields.Point(1))

    H = 0.5 * (E3 - E1) * inv(im * ω * μ)
    Z = E2 / H
    rho_a = get_appres(Z, ω)
    phi = get_phase(Z)

    return rho_a, phi
end

function SubsurfaceCore.forward(m::Tm,
        ω::T3,
        response_trans_utils::T=default_mt_tf_fns) where {
        Tm <: MTModel{<:AbstractVector, <:GridapType}, T, T3}
    if !(num_cells(m.h) == length(m.m))
        error("mesh not aligned with model values")
    end
    n = length(ω)
    ρₐ = zeros(eltype(m.m), n)
    ϕ = zeros(eltype(m.m), n)
    i = 1
    @inbounds while i <= n
        ρₐ[i], ϕ[i] = get_Z(m, ω[i])
        i += 1
    end
    f1 = response_trans_utils.ρₐ.tf
    f2 = response_trans_utils.ϕ.tf
    MTResponse{typeof(ρₐ), typeof(ϕ)}((f1.(ρₐ)), f2.(ϕ))
end
