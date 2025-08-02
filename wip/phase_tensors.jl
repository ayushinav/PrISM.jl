using Pkg

Pkg.activate("wip/.")

using LinearAlgebra
using DifferentiationInterface
using Enzyme

function R(θ_p)
    θ = θ_p *π/180
    reshape([cos(θ), -sin(θ), sin(θ), cos(θ)], 2,2)
end

x_obs = randn(2)
X_true = [0 x_obs[1]; x_obs[2] 0] 
y_obs = randn(2)
Y_true =[0 y_obs[1]; y_obs[2] 0] 

Z_true = @. X_true + im * Y_true
Φ_true = inv(X_true) * Y_true

θ_rot = rand() * 360
X_obs = R(-θ_rot) * X_true * R(θ_rot)
Y_obs = R(-θ_rot) * Y_true * R(θ_rot)

Z_obs = R(-θ_rot) * Z_true * R(θ_rot)
Z_err = 0.01 .* Z_obs
# Z_obs = @. X_obs + im * Y_obs
# diagonal elements are conjugates of each other

Φ_obs = inv(X_obs) * Y_obs
# Φ_obs = R(-θ_rot) * Φ_true * R(θ_rot)
# off diagonal elements are equal

get_β(M) = 0.5 * atan((M[1,2] - M[2,1])/(M[1,1] + M[2,2])) * 180/π


function get_β_from_pars(x)
    X = reshape(x[1:4], 2,2)
    Y = reshape(x[5:8], 2,2)

    Φ = inv(X) *Y
    get_β(Φ)
end

get_β_from_pars([X_obs[:]..., Y_obs[:]...])

pars = [X_obs[:]..., Y_obs[:]...]

gp = DifferentiationInterface.gradient(get_β_from_pars, AutoFiniteDiff(), pars)

g_pp = DifferentiationInterface.hessian(get_β_from_pars, AutoFiniteDiff(), pars)

mean_pars = [X_obs[:]..., Y_obs[:]...];
std_pars = [X_obs[:]..., Y_obs[:]...] .* 0.01;

mean_β = get_β_from_pars(mean_pars) + 0.5 * diag(g_pp)' * std_pars
std_β = sum(gp.^ 2 .* std_pars)


