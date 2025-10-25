# helper functions

import LinearAlgebra: zero

function zero(x::resp) where {resp <: AbstractGeophyResponse}
    typeof(x)([zero(getfield(x, k)) for k in fieldnames(resp)]...)
end
function zero(x::model) where {model <: AbstractGeophyModel}
    typeof(x)([zero(getfield(x, k)) for k in fieldnames(model)]...)
end

import Base: copy
function copy(x::resp) where {resp <: AbstractGeophyResponse}
    typeof(x)([copy(getfield(x, k)) for k in fieldnames(resp)]...)
end
function copy(x::model) where {model <: AbstractGeophyModel}
    typeof(x)([copy(getfield(x, k)) for k in fieldnames(model)]...)
end

# forward manipulation

# function SubsurfaceCore.forward_helper(
#         m::Type{T}, m0, vars, response_trans_utils, params) where {T <: AbstractGeophyModel}
#     model = from_nt(m, m0)
#     resp_nt = to_resp_nt(forward(model, vars, response_trans_utils))
#     return resp_nt
# end

# Only occam uses the following

function zero_abstract(m::mtresponse) where {mtresponse <: MTResponse{
        <:AbstractVector{<:Any}, <:AbstractVector{<:Any}}}
    MTResponse{AbstractVector, AbstractVector}(zero(m.ρₐ), zero(m.ϕ))
end

# function inverse(t::mtresponse; abstract=false) where {mtresponse <: MTResponse}
#     if abstract
#         return MTModel{AbstractArray{<:Any, length(size(t.ρₐ))}, # length(size(...)) => dimensionality
#             AbstractArray{<:Any, length(size(t.ρₐ))}}
#     else
#         vec_type = typeof(t.ρₐ)
#         return MTModel{vec_type, vec_type}
#     end
# end

# function forward(t::mtmodel; abstract=false) where {mtmodel <: MTModel} # {<:AbstractVector{<:Any}, <:AbstractVector{<:Any}}}
#     if abstract
#         return MTResponse{AbstractArray{<:Any, length(size(m.m))},
#             AbstractArray{<:Any, length(size(m.m))}}
#     else
#         vec_type = typeof(t.m)
#         return MTResponse{vec_type, vec_type}
#     end
# end
