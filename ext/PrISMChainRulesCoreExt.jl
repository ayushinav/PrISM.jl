module PrISMChainRulesCoreExt

using PrISM, ChainRulesCore
import ChainRulesCore: frule, rrule
import PrISM: find_c_

function ChainRulesCore.frule(
        config::ChainRulesCore.RuleConfig{>:ChainRulesCore.HasForwardsMode},
        (_, _, _, _, Δm, _), ::typeof(find_c_), f, c1, c2, ω, m::T) where {T <: RWModel}
    c = find_c_(f, c1, c2, ω, m)

    function ∇_find_c_(Δm)
        _,
        fₓ = ChainRulesCore.frule_via_ad(config, (NoTangent(), one(c)), x -> f(x, ω, m), c)
        _, fₚ = ChainRulesCore.frule_via_ad(config, (NoTangent(), Δm), m_ -> f(c, ω, m_), m)
        ∇x = -fₚ / fₓ
    end

    return c, ∇_find_c_(Δm)
end

# function ChainRulesCore.frule(config::ChainRulesCore.RuleConfig{>:ChainRulesCore.HasForwardsMode},
#     (_, _, _, _, Δm, _), 
#     ::typeof(find_c_), f, c1, c2, ω, m::T) where {T <: LWModel}
#     c = find_c_(f, c1, c2, ω, m)

#     function ∇_find_c_(Δm)
#         _, fₓ = ChainRulesCore.frule_via_ad(config, (NoTangent(), one(c)), x -> f(x, ω, m), c)
#         _, fₚ = ChainRulesCore.frule_via_ad(config, (NoTangent(), Δm),     m_ -> f(c, ω, m_), m)
#         ∇x = -fₚ / fₓ
#     end

#     return c, ∇_find_c_(Δm)
# end

end
