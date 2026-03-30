@testitem "type inference tests : MT" tags = [:inference] begin
    using BenchmarkTools, JET
    h = [100.0, 1000.0] # m
    ρ = log10.([100.0, 10.0, 1000.0]) # Ωm
    m = MTModel(ρ, h)
    T = 10 .^ (range(-3, 5; length=57))
    ω = 2π ./ T
    nω = length(T)
    @inferred MTModel(ρ, h)
    @inferred forward(m, ω)
    @test_opt forward(m, ω)
    @test_call forward(m, ω)
end

@testitem "type inference tests : Rayleigh Waves" tags = [:inference] begin
    using BenchmarkTools
    velocity_model = [
        [10.0, 7.00, 3.50, 2.00], [10.0, 6.80, 3.40, 2.00], [10.0, 7.00, 3.50, 2.00],
        [10.0, 7.60, 3.80, 2.00], [10.0, 8.40, 4.20, 2.00], [10.0, 9.00, 4.50, 2.00],
        [10.0, 9.40, 4.70, 2.00], [10.0, 9.60, 4.80, 2.00], [10.0, 9.50, 4.75, 2.00]]

    vmodel = hcat(velocity_model...)'

    h_ = vmodel[:, 1]
    vp_ = (vmodel[:, 2]) #.+ 4.5 * 1.72
    vs_ = (vmodel[:, 3]) #.+ 4.5
    density_ = vmodel[:, 4]
    m_rw = RWModel(vs_, h_[1:(end - 1)], density_, vp_)

    t = exp10.(range(0, 3, length=100))

    @inferred RWModel(vs_, h_[1:(end - 1)], density_, vp_)
    @inferred forward(m_rw, t)
    @test_opt forward(m_rw, t)
    @test_call forward(m_rw, t)

    # Performance test
    # alloc = @ballocated forward!(resp, m, ω)

end

@testitem "type inference tests : DC Resistivity" tags = [:inference] begin
    ρ = log10.([1e3, 4e3, 2e2])
    h = [100.0, 100.0]
    m = DCModel(ρ, h)
    locs = get_wenner_array(range(20, 500, length=25))

    @inferred DCModel(ρ, h)
    @inferred forward(m, locs)
    @test_opt forward(m, locs)
    @test_call forward(m, locs)
end
