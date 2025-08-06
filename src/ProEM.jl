module ProEM
using Reexport
@reexport using SubsurfaceCore
using LinearAlgebra
using LinearSolve
using NonlinearSolve
using Optimization, OptimizationOptimJL
using ProgressMeter
using Enzyme
using DifferentiationInterface
using UnPack
using InteractiveUtils
using Distributions
using Turing
using Gridap
using GridapGmsh
using GridapMakie
using Gridap.Fields
using GLMakie
import Base: show

import SubsurfaceCore: forward, forward_helper

include("models/mt/types.jl")
include("models/mt/forward.jl")
include("models/mt/forward_fem.jl")
include("models/pretty_printing.jl")

include("plots/fem_1d.jl")
include("utils.jl")
include("inverse/utils.jl")
include("inverse/occam.jl")
include("inverse/inv.jl")
include("inverse/nl_inv.jl")
include("inverse/opt_inv.jl")
include("probabilistic/models/mt.jl")
include("probabilistic/utils.jl")
include("probabilistic/rto.jl")

include("../extension/fem_utils.jl")
include("../extension/fem1d_gmsh_generator.jl")

# geophysics
export MTModel, MTResponse

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

end
