using Pkg
Pkg.activate("wip/.")

using JLD2
using Enzyme
using DifferentiationInterface
using ForwardDiff
using FiniteDiff
using LinearAlgebra

z_arr = jldopen("../../kivu/tfs_rotated/unrotated.jld2")["z_arr"]
zerr_arr = jldopen("../../kivu/tfs_rotated/unrotated.jld2")["zerr_arr"]
f_arr = jldopen("../../kivu/tfs_rotated/unrotated.jld2")["f_arr"]

get_β(M) = 0.5 * atan((M[1,2] - M[2,1])/(M[1,1] + M[2,2])) * 180/π
function R(θ_p)
    θ = θ_p *π/180
    reshape([cos(θ), -sin(θ), sin(θ), cos(θ)], 2,2)
end

function get_β_from_pars(x)
    X = reshape(x[1:4], 2,2)
    Y = reshape(x[5:8], 2,2)

    Φ = inv(X) *Y
    get_β(Φ)
end


err_Φ_ij(μx, μy, σx², σy²) = sqrt((μx^2* σy² + μy^2* σx²)/ (μx^4))
err_Φ_ij(μx, μy, σ²)= err_Φ_ij(μx, μy, σ², σ²)


function get_β_err(Φ, ΣΦ)
    k = (Φ[1,2]- Φ[2,1])/(Φ[1,1]+ Φ[2,2])
    σ = 0.5 * 180/π * sqrt(
        1/(1+k^2)^2 * inv(Φ[1,1]+ Φ[2,2])^4 * (
            (Φ[1,1] + Φ[2,2])^2 * (ΣΦ[1,2] + ΣΦ[2,1]) +
            (Φ[1,2] + Φ[2,1])^2 * (ΣΦ[1,1] + ΣΦ[2,2])
        )
    )

    return σ
    
end


function get_θ_from_pars(P_sym_)
    # X = reshape(x[1:4], 2,2)
    # Y = reshape(x[5:8], 2,2)

    # P = inv(X) *Y
    # β = get_β(P)

    # P_sym = R(-β) * P * R(-β)


    P_sym = reshape(P_sym_, 2, 2)
    l, q = eigen(P_sym)

    v1 = q[:, 1]
    # v2 = q[:, 2]
    θ = atan(v1[2]/v1[1]) * 180/π
    return θ
end

x_arr = real.(z_arr)
y_arr = imag.(z_arr)

Φ_arr = [zero(x_arr[i]) for i in eachindex(z_arr)]
β_arr = [zero(x_arr[i][:,1,1]) for i in eachindex(z_arr)]
βerr_arr = [zero(x_arr[i][:,1,1]) for i in eachindex(z_arr)]
Φsym_arr = [zero(x_arr[i]) for i in eachindex(z_arr)]
eig_arr = [zero(x_arr[i]) for i in eachindex(z_arr)]

isym_arr = [zero(x_arr[i][:,1,1]) for i in eachindex(z_arr)]
eig_sym_arr = [zero(x_arr[i][:,1,1]) for i in eachindex(z_arr)]

θ_arr = [zero(x_arr[i][:,1,1]) for i in eachindex(z_arr)]
θerr_arr = [zero(x_arr[i][:,1,1]) for i in eachindex(z_arr)]

new_z_arr = [zero(z_arr[i]) for i in eachindex(z_arr)]


for i in eachindex(z_arr)
    for j in axes(z_arr[i], 1)
        P = inv(x_arr[i][j,:,:]) * y_arr[i][j,:,:]
        Φ_arr[i][j,:,:] .= P
        mean_pars = [x_arr[i][j,:,:][:]..., y_arr[i][j,:,:]...]
        std_pars = [zerr_arr[i][j,:,:]..., zerr_arr[i][j,:,:]...]
        
        P_err = err_Φ_ij.(x_arr[i][j,:,:], y_arr[i][j,:,:], zerr_arr[i][j,:,:].^2)

        β = get_β(P)
        β_arr[i][j] = β
        βerr = get_β_err(P, P_err .^ 2)
        βerr_arr[i][j] = βerr

        P_sym = R(-β) * P * R(-β)

        K = kron(R(β), R(-β))
        ΣP_vec = diagm(P_sym[:]).^2
        ΣP_sym_vec = K * ΣP_vec * K'
        
        # Φsym_arr[i][j,:,:] .= P_sym
        l, q = eigen(P_sym)

        # v1 = q[:, 1]
        # v2 = q[:, 2]

        gp = DifferentiationInterface.gradient(get_θ_from_pars, AutoFiniteDiff(), P_sym[:])
        g_pp = DifferentiationInterface.hessian(get_θ_from_pars, AutoFiniteDiff(), P_sym[:])

        # θ = atan(v1[2]/v1[1]) * 180/π
        # θ_arr[i][j] = θ

        θ = get_θ_from_pars(P_sym[:]) #+ 0.5 * diag(g_pp)' * std_pars
        θerr =  sqrt(sum(gp.^ 2 .* diag(ΣP_sym_vec)))
        θ_arr[i][j] = θ
        θerr_arr[i][j] = θerr
        
        eig_sym_arr[i][j] = norm(R(θ) * P_sym * R(-θ) .- diagm(l), 1)

    end
end

[norm(eig_sym_arr[i], 1) for i in eachindex(z_arr)] .< 1e-10




