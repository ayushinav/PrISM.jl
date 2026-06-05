# Automatic differentiation

```@setup ad_compats
using DifferentiationInterface
using Enzyme
using PrISM
using ForwardDiff
using PrettyTables
using BenchmarkTools
using Printf
```

## Optimization

AD compatibility for optimization can be summarized in the following table. This affects deterministic inversion and RTO-TKO. We also show the runtimes for a 50-layered model in all the cases when jacobians are correct :

```@raw html
<details closed><summary>Code for the tables</summary>
```

```@example ad_compats
model_types = [MTModel, RWModel, LWModel]

sigmoid_lw(x) = SubsurfaceCore.sigmoid(x, 3.0, 5.0)
models = ((; m=randn(50) .* 1.0 .+ 3, h=fill(100.0, 49)),
    (; m=rand(50) .* 20e-1 .+ 3.0, h=fill(100.0, 49), vp=fill(7.5, 50), ρ=fill(2.5, 50)),
    (; m=sigmoid_lw.(cumsum(rand(50))), h=fill(100.0, 49), ρ=fill(2.5, 50)),
    (; m=randn(50) .* 0.01 .+ 2, h=fill(100.0, 49)))

dc_locs = get_wenner_array(range(20, 500; length=25))
vars = (10.0 .^ collect(-3:0.1:1), 10.0 .^ collect(-1:0.1:1),
    10.0 .^ collect(-1:0.1:1), dc_locs)

ADTYPES = (AutoFiniteDiff(),)
Cache_DI = DifferentiationInterface.Cache
Constant_DI = DifferentiationInterface.Constant
jac_baseline = []

for ik in eachindex(model_types), adtype_ in ADTYPES

    model_ref = from_nt(model_types[ik], models[ik])
    model_cache = deepcopy(model_ref)
    m = model_ref.m
    vars_ = vars[ik]

    resp = forward(model_ref, vars_)
    resp_cache_ = deepcopy(resp)
    response_fields_ = propertynames(resp)

    kk = ntuple(i -> no_tf, length(propertynames(resp)))
    response_trans_utils_ = (; zip(propertynames(resp), kk)...)

    model_trans_utils_ = no_tf
    model_type_ = model_types[ik]
    @show model_type_

    rvec_ = vcat(values(deepcopy(to_nt(resp)))...)
    params = default_params(model_type_)

    prep_j = prepare_jacobian(
        PrISM.wrapper_DI!, rvec_, adtype_, m, Cache_DI(model_cache), Constant_DI(model_ref),
        Cache_DI(resp_cache_), Constant_DI(vars_), Constant_DI(response_fields_),
        Constant_DI(model_type_), Constant_DI(model_trans_utils_),
        Constant_DI(response_trans_utils_), Constant_DI(params))

    jacobian_ = zeros(length(rvec_), length(m));

    DifferentiationInterface.jacobian!(
        PrISM.wrapper_DI!, rvec_, jacobian_, prep_j, adtype_, m,
        Cache_DI(model_cache), Constant_DI(model_ref), Cache_DI(resp_cache_),
        Constant_DI(vars_), Constant_DI(response_fields_),
        Constant_DI(model_type_), Constant_DI(model_trans_utils_),
        Constant_DI(response_trans_utils_), Constant_DI(params))

    push!(jac_baseline, jacobian_)
end

ADTYPES = (
    AutoFiniteDiff(), AutoForwardDiff(), AutoEnzyme(; mode=set_runtime_activity(Reverse)),
    AutoEnzyme(; mode=set_runtime_activity(Forward)))
#

flag_works = fill(true, (length(model_types), length(ADTYPES)))
flag_acc = fill(false, (length(model_types), length(ADTYPES)))
bm_times = fill("✘", (length(model_types), length(ADTYPES)))

for ik in eachindex(model_types), i_ad in eachindex(ADTYPES)

    adtype_ = ADTYPES[i_ad]

    model_ref = from_nt(model_types[ik], models[ik])
    model_cache = deepcopy(model_ref)
    m = model_ref.m
    vars_ = vars[ik]

    resp = forward(model_ref, vars_)
    resp_cache_ = deepcopy(resp)
    response_fields_ = propertynames(resp)

    kk = ntuple(i -> no_tf, length(propertynames(resp)))
    response_trans_utils_ = (; zip(propertynames(resp), kk)...)

    model_trans_utils_ = no_tf

    model_type_ = model_types[ik]

    rvec_ = vcat(values(deepcopy(to_nt(resp)))...)
    params = default_params(model_type_)

    try
        prep_j = prepare_jacobian(
            PrISM.wrapper_DI!, rvec_, adtype_, m, Cache_DI(model_cache),
            Constant_DI(model_ref), Cache_DI(resp_cache_),
            Constant_DI(vars_), Constant_DI(response_fields_),
            Constant_DI(model_type_), Constant_DI(model_trans_utils_),
            Constant_DI(response_trans_utils_), Constant_DI(params))

        jacobian_ = zeros(length(rvec_), length(m))

        DifferentiationInterface.jacobian!(PrISM.wrapper_DI!, rvec_, jacobian_, adtype_, m,
            Cache_DI(model_cache), Constant_DI(model_ref), Cache_DI(resp_cache_),
            Constant_DI(vars_), Constant_DI(response_fields_),
            Constant_DI(model_type_), Constant_DI(model_trans_utils_),
            Constant_DI(response_trans_utils_), Constant_DI(params))

        flag_acc[ik, i_ad] = isapprox(jacobian_, jac_baseline[ik]; rtol=0.1)

        time_ = @belapsed begin
            DifferentiationInterface.jacobian!(
                PrISM.wrapper_DI!, $rvec_, $jacobian_, $adtype_, $m,
                Cache_DI($model_cache), Constant_DI($model_ref), Cache_DI($resp_cache_),
                Constant_DI($vars_), Constant_DI($response_fields_),
                Constant_DI($model_type_), Constant_DI($model_trans_utils_),
                Constant_DI($response_trans_utils_), Constant_DI($params))
        end
        bm_times[ik, i_ad] = @sprintf("%.3f ms", time_*1e3)

    catch
    end
end

adtypes_string = split("FiniteDiff ForwardDiff Enzyme:Reverse Enzyme:Forward", " ")
bm_times[.!flag_acc] .= "✔"
nothing # hide
```

```@raw html
</details>
```

```@example ad_compats
pretty_table(bm_times; row_labels=[string.(model_types)...], header=adtypes_string) # hide
```

## Turing AD compatibility

Compatibility with using AD inside samplers such as Hamiltonian MCMC and No U-Turn (NUTS) can be summarized in the following table. We also report runtimes for 50 layered model whenever the models execute without errors.

```@raw html
<details closed><summary>Code for the table</summary>
```

```@example ad_compats
using Distributions, Turing, Printf
using DynamicPPL.TestUtils.AD: run_ad

modelD_types_t = [
    MTModelDistribution, RWModelDistribution, LWModelDistribution, DCModelDistribution]
respD_t = [MTResponseDistribution(normal_dist, normal_dist),
    SurfaceWaveResponseDistribution(normal_dist),
    SurfaceWaveResponseDistribution(normal_dist), DCResponseDistribution(normal_dist)]

sigmoid_lw(x) = SubsurfaceCore.sigmoid(x, 3.0, 5.0)
true_models_t = ((; m=randn(50) .* 1.0 .+ 3, h=fill(100.0, 49)),
    (; m=rand(50) .* 20e-1 .+ 3.0, h=fill(100.0, 49), vp=fill(7.5, 50), ρ=fill(2.5, 50)),
    (; m=sigmoid_lw.(cumsum(rand(50))), h=fill(100.0, 49), ρ=fill(2.5, 50)),
    (; m=randn(50) .* 0.01 .+ 2, h=fill(100.0, 49)))
vars_t = [10.0 .^ collect(-3:0.1:1), 10.0 .^ collect(0:0.1:2), 10.0 .^ collect(-1:0.1:1)]

ADTYPES_turing = (AutoForwardDiff(), AutoEnzyme(; mode=set_runtime_activity(Reverse)))
adtypes_turing_string = ["ForwardDiff", "Enzyme:Reverse"]

turing_results = fill("✘", (length(model_types), length(ADTYPES_turing)))

for ik in eachindex(model_types), i_ad in eachindex(ADTYPES_turing)

    adtype_ = ADTYPES_turing[i_ad]

    model = from_nt(model_types[ik], true_models_t[ik])
    resp = forward(model, vars_t[ik])

    err_resp = copy(resp)
    for k in propertynames(resp)
        setproperty!(err_resp, k, getproperty(err_resp, k) .* 0.05)
    end

    mdist = from_nt(modelD_types_t[ik],
        (; true_models_t[ik]...,
            m=product_distribution([Uniform(mi * 0.8, mi * 1.2) for mi in model.m])))

    m_cache = mcmc_cache(mdist, respD_t[ik])
    dppl_object, _ = get_stochastic_inverse_model(resp, err_resp, vars_t[ik], m_cache)

    try
        run_ad(dppl_object, adtype_; atol=1e-6, params=rand(mdist.m))
    catch
    end

    turing_results[ik, i_ad] = try
        m_ = rand(mdist.m)
        elapsed = @belapsed run_ad($dppl_object, $adtype_; atol=1e-6, params=($m_))
        @sprintf("%.2f ms", elapsed*1e3)
    catch
        "✘"
    end
end
```

```@raw html
</details>
```

```@example ad_compats
pretty_table(turing_results; row_labels=[string.(model_types)...], header=adtypes_turing_string) # hide
```

!!! compat

    Gradients and jacobians for `LWModels` are correct only when the shear wave velocities increase with depth.

### Reproducibility

```@example ad_compats
println("The above benchmarks were obtained on $(Sys.cpu_info()[1].model)") # hide
# println(versionfo()) # hide
```
