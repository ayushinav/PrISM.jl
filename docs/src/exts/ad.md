# Automatic differentiation

```@setup ad_compats
using DifferentiationInterface
using Enzyme
using PrISM
using ForwardDiff
using PrettyTables
```

The current functionality of automatic differentiation (AD) can be summarized in the following tables :

```@raw html
<details closed><summary>Code for the tables</summary>
```

```@example ad_compats
model_types = [MTModel, RWModel, LWModel] #, DCModel]
modelD_types = [MTModelDistribution, RWModelDistribution, LWModelDistribution]
respD = [MTResponseDistribution(normal_dist, normal_dist),
    SurfaceWaveResponseDistribution(normal_dist),
    SurfaceWaveResponseDistribution(normal_dist)]

models = ((; m=randn(50) .* 0.01 .+ 2, h=fill(100.0, 49)),
    (; m=randn(50) .* 10e-3 .+ 3.5, h=fill(100.0, 49), vp=fill(7.5, 50), ρ=fill(2.5, 50)),
    (; m=randn(50) .* 10e-3 .+ 3.5, h=fill(100.0, 49), ρ=fill(2.5, 50)),
    (; m=randn(50) .* 0.01 .+ 2, h=fill(100.0, 49)))
dc_locs = get_wenner_array(range(20, 500; length=25))
vars = (
    10.0 .^ collect(-3:0.1:1), 10.0 .^ collect(0:0.1:3), 10.0 .^ collect(0:0.1:3), dc_locs)

# samplers = [NUTS, SliceSampler]
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

    # kk = ntuple(i -> no_tf, length(propertynames(model_ref)))
    model_trans_utils_ = no_tf #(; zip(propertynames(model_ref), kk )...)

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

    # @show respon

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

for ik in eachindex(model_types), i_ad in eachindex(ADTYPES)
    # ik, i_ad = 3,2
    # @show model_types[ik]
    # @show ADTYPES[i_ad]

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

    # kk = ntuple(i -> no_tf, length(propertynames(model_ref)))
    model_trans_utils_ = no_tf #(; zip(propertynames(model_ref), kk )...)

    model_type_ = model_types[ik]

    rvec_ = vcat(values(deepcopy(to_nt(resp)))...)
    params = default_params(model_type_)

    # flag_1 = true
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

        flag_acc[ik, i_ad] = all(isapprox.(jacobian_, jac_baseline[ik]; atol=0.4))

        # @show extrema(jac_baseline[ik])
        # @show extrema(jacobian_)
        # @show "============"

    catch
        flag_works[ik, i_ad] = false
    end

    # @show flag_works[ik, i_ad]
    # @show flag_acc[ik, i_ad]

    # break
end

adtypes_string = split("FiniteDiff ForwardDiff Enzyme:Reverse Enzyme:Forward", " ")
function get_Table(mask)
    tick_table = fill("✘", size(mask))
    tick_table[mask] .= "✔"
    pretty_table(tick_table; row_labels=[string.(model_types)...], header=adtypes_string)
end
```

```@raw html
</details>
```

## Stability

```@example ad_compats
work_table = get_Table(flag_works) # hide
```

## Accuracy

```@example ad_compats
acc_table = get_Table(flag_acc) # hide
```
