"""
    get_schlumberger_array(a, n)
returns electrode spacings for Schlumberger array :

M ←na→ A ←a→ B ←na→ N

## Arguments
 - `a` : electrode spacing AB
 - `n` : vector of numbers corresponding to spacing

## Usage

"""
function get_schlumberger_array(a, n::T) where T <: AbstractVector
    srcs = hcat([[-i*a - a/2, i*a + a/2] for i in n]...)'
    recs = [-a/2, a/2]
    locs = (; recs, srcs)
end

"""
    get_wenner_array(a)
returns electrode spacings for Wenner array :

M ←a→ A ←a→ B ←a→ N

## Arguments
 - `a` : vector of electrode spacing AB

## Usage

TODO
"""
function get_wenner_array(a::T) where T <: AbstractVector
    recs = hcat(-a./2, a./2)
    srcs = hcat(-a.* 1.5f0, a.* 1.5f0)
    locs = (; recs, srcs)
end

