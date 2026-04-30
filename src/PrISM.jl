module PrISM
using Reexport
@reexport using SubsurfaceCore
using LinearAlgebra
using StaticArrays
using LinearSolve
using ProgressMeter
using DifferentiationInterface
using UnPack # TODO
using InteractiveUtils # TODO 
using PrettyTables
using DataInterpolations
using Distributions
import Base: show

# import DifferentiationInterface: recursive_similar

import SubsurfaceCore: forward, forward_helper, get_scales, get_labels

include("models/filters.jl")

include("models/mt/types.jl")
include("models/mt/utils.jl")
include("models/mt/forward.jl")

include("models/surface_waves/types.jl")
include("models/surface_waves/utils.jl")
include("models/surface_waves/forward.jl")

include("models/dc/types.jl")
include("models/dc/utils.jl")
include("models/dc/forward.jl")
include("models/dc/arrays.jl")

include("models/pretty_printing.jl")

include("utils.jl")
include("inverse/utils.jl")
include("inverse/occam.jl")
include("inverse/inv.jl")
include("probabilistic/models/mt.jl")
include("probabilistic/models/surface_waves.jl")
include("probabilistic/models/dc.jl")
include("probabilistic/utils.jl")
include("probabilistic/rto.jl")
include("plots/utils.jl")

# geophysics
export MTModel, MTResponse
export RWModel, LWModel, SurfaceWaveResponse
export DCModel, DCResponse
export get_schlumberger_array, get_wenner_array

# forward
export get_Z, get_appres, get_phase, forward!, forward

#inverse
export occam_cache, Occam, nl_cache, NonlinearAlg, opt_cache, OptAlg, linsolve!, occam_step!
export inverse!
export ∂, χ², linear_utils, inverse_utils

# probabilistic
export rto_cache

## geophysics
export MTModelDistribution, MTResponseDistribution
export RWModelDistribution, SurfaceWaveResponseDistribution
export LWModelDistribution
export DCModelDistribution, DCResponseDistribution

# TODO
function SubsurfaceCore.forward_helper(
        m::Type{T}, m0, vars, response_trans_utils, params) where {T <: AbstractGeophyModel}
    model = from_nt(m, m0)
    resp_nt = to_resp_nt(forward(model, vars, params))
    for k in propertynames(resp_nt)
        broadcast!(response_trans_utils[k].tf, resp_nt[k], resp_nt[k])
    end
    return resp_nt
end
end
