module PrISMForwardDiffExt

using PrISM, ChainRulesCore, ForwardDiff, ForwardDiffChainRules
import PrISM: find_c_

@ForwardDiff_frule find_c_(f, c1, c2, m::AbstractVector{<:ForwardDiff.Dual}, ω)

end
