module ProEM
using Reexport
@reexport using SubsurfaceCore
using LinearAlgebra
using StaticArrays
using LinearSolve
using NonlinearSolve
using Optimization, OptimizationOptimJL
using ProgressMeter
using Enzyme
using DifferentiationInterface
using UnPack
using InteractiveUtils
using DataInterpolations
using Distributions
using Turing
import Base: show

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

include("models/pretty_printing.jl")

include("utils.jl")
include("inverse/utils.jl")
include("inverse/occam.jl")
include("inverse/inv.jl")
include("inverse/nl_inv.jl")
include("inverse/opt_inv.jl")
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

end
