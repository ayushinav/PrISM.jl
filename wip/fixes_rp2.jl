using Pkg
Pkg.activate(".")

using MT
using Test
using Distributions
using Turing
using BenchmarkTools

n1 = multi_phase_modelType(SEO3, Sifre2014, Zhang2012, HS_minus_multi_phase())

T = [1200f0] #collect(1000:100:1400) .+ 273.0f0
Ch2o_ol = 1.0f0
Ch2o_m = 1000.0f0
Cco2_m = 1000.0f0
Ch2o_opx = 100.0f0
Cco2_m = 10.0f0
ϕ = [0.1f0, 0.2f0]

ps_nt = (; ϕ=ϕ, T=T, Ch2o_ol=Ch2o_ol, Ch2o_m=Ch2o_m, Cco2_m=Cco2_m, Ch2o_opx)

m1 = n1(ps_nt)

@benchmark forward(m1, [])
#= 

BenchmarkTools.Trial: 10000 samples with 10 evaluations per sample.
 Range (min … max):  1.333 μs … 697.312 μs  ┊ GC (min … max): 0.00% … 97.66%
 Time  (median):     1.392 μs               ┊ GC (median):    0.00%
 Time  (mean ± σ):   1.636 μs ±   9.003 μs  ┊ GC (mean ± σ):  4.16% ±  0.98%

   █▅                                                          
  ▄██▆▃▃▃▂▂▂▂▁▂▂▂▂▂▁▂▂▂▁▂▂▂▂▂▂▂▂▁▂▂▂▂▂▂▂▂▂▂▂▂▂▂▂▂▂▂▂▂▂▂▂▂▂▂▁▂ ▂
  1.33 μs         Histogram: frequency by time        3.03 μs <

 Memory estimate: 2.03 KiB, allocs estimate: 38.

=#

forward(m1, [])

@inferred forward(m1, [])

@report_opt forward(m1, [])
@report_call forward(m1, [])

@benchmark forward(m1, [])

#=
BenchmarkTools.Trial: 10000 samples with 96 evaluations per sample.
 Range (min … max):  788.188 ns … 93.708 μs  ┊ GC (min … max):  0.00% … 97.99%
 Time  (median):     829.865 ns              ┊ GC (median):     0.00%
 Time  (mean ± σ):   978.990 ns ±  3.034 μs  ┊ GC (mean ± σ):  11.77% ±  3.77%

  ▄██▇▇▆▆▅▄▃▂▁▁▁▁▁                                             ▂
  ██████████████████▇██▇▇▇▆▆▆▆▅▇▇▆▆▇▅▅▅▂▆▄▄▅▅▅▅▅▄▅▄▄▃▄▃▂▄▃▄▂▄▅ █
  788 ns        Histogram: log(frequency) by time      1.51 μs <

 Memory estimate: 1.70 KiB, allocs estimate: 31.
=#


function return_broadcasting_function(phi, mix, ::Val{2})
    function fn(m1, m2)
        mix_models((m1, m2), phi, mix)
    end
end

function return_broadcasting_function(phi, mix, ::Val{3})
    function fn(m1, m2, m3)
        mix_models((m1, m2, m3), phi, mix)
    end
end

function return_broadcasting_function(phi, mix, ::Val{4})
    function fn(m1, m2, m3, m4)
        mix_models((m1, m2, m3, m4), phi, mix)
    end
end

function return_broadcasting_function(phi, mix, ::Val{5})
    function fn(m1, m2, m3, m4, m5)
        mix_models((m1, m2, m3, m4, m5), phi, mix)
    end
end

function broadcast_helper_(m_tup, phi, mix, ::Val{2})
    broadcast((m1, m2) -> mix_models((m1, m2), model.ϕ, model.mix), m_tup[1], m_tup[2])
end

function broadcast_helper_(m_tup, phi, mix, ::Val{3})
    broadcast((m1, m2, m3) -> mix_models((m1, m2, m3), model.ϕ, model.mix), m_tup[1], m_tup[2], m_tup[3])
end

function broadcast_helper_(m_tup, phi, mix, ::Val{4})
    broadcast((m1, m2, m3, m4) -> mix_models((m1, m2, m3, m4), model.ϕ, model.mix), m_tup[1], m_tup[2], m_tup[3], m_tup[4])
end

function broadcast_helper_(m_tup, phi, mix, ::Val{5})
    broadcast((m1, m2, m3, m4, m5) -> mix_models((m1, m2, m3, m4, m5), model.ϕ, model.mix), m_tup[1], m_tup[2], m_tup[3], m_tup[4], m_tup[5])
end


arr_ = rand(10)
fnames_ = rand(Bool, 10)

function filter_type_stable(arr, fnames)
    arr[fnames]
end

@inferred filter_type_stable(arr_, fnames_)

400e-9 * 1e9 * 100 / 3600