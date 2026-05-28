using PrISM, DifferentiationInterface, Enzyme, LinearAlgebra, Test
using DifferentiationInterface: AutoEnzyme
Cache_DI = DifferentiationInterface.Cache
Constant_DI = DifferentiationInterface.Constant

model_types = [MTModel, RWModel, LWModel]

true_models = ((; m=randn(50) .* 1.0 .+ 2, h=fill(100.0, 49)),
    (; m=rand(50) .* 20e-1 .+ 3.0, h=fill(100.0, 49), vp=fill(7.5, 50), ρ=fill(2.5, 50)),
    (; m=rand(50) .* 20e-1 .+ 3, h=fill(100.0, 49), ρ=fill(2.5, 50)),
    (; m=randn(50) .* 0.01 .+ 2, h=fill(100.0, 49)))

vars = [10.0 .^ collect(-3:0.1:1), 10.0 .^ collect(0:0.1:2), 10.0 .^ collect(-1:0.1:1)]

adtype_baseline = AutoFiniteDiff()
ADTYPES = (AutoEnzyme(; mode=Enzyme.set_runtime_activity(Reverse)),
    AutoEnzyme(; mode=Enzyme.set_runtime_activity(Forward)))

@testset "$(model_types[ik]) : $(adtype.mode)" for ik in eachindex(model_types),
    adtype in ADTYPES

    model_ref = from_nt(model_types[ik], true_models[ik])
    vars_ = vars[ik]

    resp = PrISM.forward(model_ref, vars_)
    response_fields_ = propertynames(resp)
    kk = ntuple(i -> no_tf, length(propertynames(resp)))
    response_trans_utils_ = (; zip(propertynames(resp), kk)...)
    model_type_ = model_types[ik]
    params = default_params(model_type_)

    tup_DI = (
        Cache_DI(deepcopy(model_ref)), Constant_DI(model_ref), Cache_DI(deepcopy(resp)),
        Constant_DI(vars_), Constant_DI(response_fields_), Constant_DI(model_type_),
        Constant_DI(no_tf), Constant_DI(response_trans_utils_), Constant_DI(params))

    m = deepcopy(model_ref.m)
    rvec_ = vcat(values(deepcopy(to_nt(resp)))...)

    prep_baseline = prepare_jacobian(
        PrISM.wrapper_DI!, rvec_, adtype_baseline, m, tup_DI...)
    jac_baseline = zeros(length(rvec_), length(m))
    jacobian!(PrISM.wrapper_DI!, rvec_, jac_baseline,
        prep_baseline, adtype_baseline, m, tup_DI...)

    prep_j = prepare_jacobian(PrISM.wrapper_DI!, rvec_, adtype, m, tup_DI...)
    jac = zeros(length(rvec_), length(m))
    jacobian!(PrISM.wrapper_DI!, rvec_, jac, prep_j, adtype, m, tup_DI...)

    @test isapprox(jac, jac_baseline; atol=0.1)
end
