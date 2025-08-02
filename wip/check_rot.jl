using LinearAlgebra

msk = reshape([0,1,1,0],2,2)
# X_obs = msk .* rand(2,2) * 100 
# Y_obs = msk .* rand(2,2) * 100


function R(θ_p)
    θ = θ_p *π/180
    reshape([cos(θ), -sin(θ), sin(θ), cos(θ)], 2,2)
end

##  ====== 1D =========
x_obs = rand()
X_obs = [0 x_obs; -x_obs 0] 
y_obs = rand()
Y_obs =[0 y_obs; -y_obs 0] 

Z_obs = @. X_obs + im * Y_obs

Φ_obs = inv(X_obs) * Y_obs
# equals scale identity matrix

get_β(M) = 0.5 * atan((M[1,2] - M[2,1])/(M[1,1] + M[2,2])) * 180/π
β = get_β(Φ_obs)
# β = 0

Φ_sym = R(-β) * Φ_obs * R(-β)
# Φ_sym is the same as Φ_obs
u, k, vt = svd(Φ_sym)
#=
u, vt are identity matrices, diagm(k) is equal to Φ_obs
=#


##  ====== 2D =========
x_obs = rand(2)
X_obs = [0 x_obs[1]; -x_obs[2] 0] 
y_obs = rand(2)
Y_obs =[0 y_obs[1]; -y_obs[2] 0] 

Z_obs = @. X_obs + im * Y_obs

Φ_obs = inv(X_obs) * Y_obs
# Z_21 and -Z_21 yield the same Φ_obs

get_β(M) = 0.5 * atan((M[1,2] - M[2,1])/(M[1,1] + M[2,2])) * 180/π
β = get_β(Φ_obs)
# β = 0 

Φ_sym = R(-β) * Φ_obs * R(-β)
# Φ_sym is the same as Φ_obs
u, k, vt = svd(Φ_sym)
#=
u, vt are identity matrices, diagm(k) is equal to Φ_obs
=#


##  ====== 3D =========
x_obs1 = rand(2)
x_obs2 = rand(2) .* 0.1
X_obs = [x_obs2[1] x_obs1[1]; -x_obs1[2] x_obs2[2]] 
y_obs1 = rand(2)
y_obs2 = rand(2) .* 0.1
Y_obs = [y_obs2[1] y_obs1[1]; -y_obs1[2] y_obs2[2]] 

Z_obs = @. X_obs + im * Y_obs

Φ_obs = inv(X_obs) * Y_obs
# Z_21 and -Z_21 yield the same Φ_obs

get_β(M) = 0.5 * atan((M[1,2] - M[2,1])/(M[1,1] + M[2,2])) * 180/π
β = get_β(Φ_obs)
# β = 0 

Φ_sym = R(-β) * Φ_obs * R(-β)
# Φ_sym is the same as Φ_obs
u, k, vt = svd(Φ_sym)
#=
u, vt are identity matrices, diagm(k) is equal to Φ_obs
=#

################################ Adding rotation ################################

##  ====== 1D =========
x_obs = rand()
X_true = [0 x_obs; -x_obs 0] 
y_obs = rand()
Y_true =[0 y_obs; -y_obs 0] 

Z_true = @. X_true + im * Y_true
Φ_true = inv(X_true) * Y_true

θ_rot = rand() * 360
X_obs = R(-θ_rot) * X_true * R(θ_rot)
Y_obs = R(-θ_rot) * Y_true * R(θ_rot)
Φ1 = inv(X1) * Y1

Z_obs = R(-θ_rot) * Z_true * R(θ_rot)
# Z_obs = @. X_obs + im * Y_obs
# still an antidiagonal matrix

Φ_obs = inv(X_obs) * Y_obs
# Φ_obs = R(-θ_rot) * Φ_true * R(θ_rot)
# scaled identity matrix

get_β(M) = 0.5 * atan((M[1,2] - M[2,1])/(M[1,1] + M[2,2])) * 180/π
β = get_β(Φ_obs)
# β = 0

Φ_sym = R(-β) * Φ_obs * R(-β)
# Φ_sym is the same as Φ_obs 

u, k, vt = svd(Φ_sym)
#=
u, vt are identity matrices, diagm(k) is equal to Φ_obs = Φ_sym
=#


##  ====== 2D =========
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
# Z_obs = @. X_obs + im * Y_obs
# diagonal elements are conjugates of each other

Φ_obs = inv(X_obs) * Y_obs
# Φ_obs = R(-θ_rot) * Φ_true * R(θ_rot)
# off diagonal elements are equal

get_β(M) = 0.5 * atan((M[1,2] - M[2,1])/(M[1,1] + M[2,2])) * 180/π
β = get_β(Φ_obs)
# β = 0

Φ_sym = R(-β) * Φ_obs * R(-β)
# Φ_sym is the same as Φ_obs 

λ, q = eigen(Φ_sym)

v1 = q[:, 1]
v2 = q[:, 2]

θ_rec2 = angle(v1[1] - im * v1[2]) * 180/π
θ_rec1 = angle(-v1[1] + im * v1[2]) * 180/π

@show norm(R(180 - θ_rec2) * Φ_sym * R(180 - -θ_rec2) .- diagm(λ), 1)
@show norm(R(180 - θ_rec1) * Φ_sym * R(180 - -θ_rec1) .- diagm(λ), 1)

θ_rec2 = atan(v1[2]/v1[1]) * 180/π

@show norm(R(θ_rec2) * Φ_sym * R(-θ_rec2) .- diagm(λ), 1)


if norm(R(θ_rec2) * Φ_sym * R(-θ_rec2) .- diagm(λ), 1) ≈ 0.
    θ_rec2 = angle(v1[1] - im * v1[2]) * 180/π
else
    θ_rec2 = 180 - angle(v1[1] - im * v1[2]) * 180/π
    @show norm(R(θ_rec2) * Φ_sym * R(-θ_rec2) .- diagm(λ), 1)
end


Φ_true
Φ_rec = R(θ_rec2) * Φ_sym * R(-θ_rec2) ## Φ_true

Φ_pred = R(β) * R(-θ_rec2) * Φ_rec * R(θ_rec2) * R(β) ## Φ_obs
Φ_obs

Z_pred = R(θ_rec2) * Z_obs * R(-θ_rec2) ## Z_true
Z_true

# if v1[1] == v2[2]
#     @show v1[2] == -v2[1]
#     θ_rec2 = angle(v1[1] - im * v1[2]) * 180/π
#     θ_rec1 = angle(-v1[1] + im * v1[2]) * 180/π


# else
#     @show v1[2] == v2[1]
#     @show "reverse"
#     θ_rec2 = angle(-v1[1] + im * v1[2]) * 180/π
# end

θ_rot
θ_rec2

Φ_true
Φ_rec = R(θ_rec2) * Φ_sym * R(-θ_rec2) ## Φ_true

Φ_pred = R(β) * R(-θ_rec2) * Φ_rec * R(θ_rec2) * R(β) ## Φ_obs
Φ_obs

Z_pred = R(θ_rec2) * Z_obs * R(-θ_rec2) ## Z_true
# Z_pred = R(θ_rec+ 180) * Z_obs * R(-θ_rec2- 180) ## Z_true
Z_true

atan.(imag.(Z_pred) ./ real.(Z_pred)).* 180/π
atan.(imag.(Z_true) ./ real.(Z_true)).* 180/π


R(180- θ_rec2) * Z_obs * R(180 - -θ_rec2)

@show θ_rot
Z_pred .- Z_true



#=
diagm(k) is equal to Φ_obs = Φ_sym
u = vt
=#

θ_grid = (0.:0.001:359.);
errs = zero(θ_grid);

for i in eachindex(θ_grid)
    c = [cos((θ_grid[i] + β/2)*π/180); sin((θ_grid[i] + β/2)*π/180)]
    p = Φ_obs *c
    c2 = [cos((θ_grid[i] - β/2)*π/180); sin((θ_grid[i] - β/2)*π/180)]
    # errs[i] = norm(p./c2)
    errs[i] = (p[1]/(c2[1] + 1e-4) - p[2]/(c2[2] + 1e-4))^2
end

argmin(errs)
θ_rec = θ_grid[argmin(errs)]




#= SIMPLER (?), but need to check in 3D case

c(ω) || Φ_sym *c(ω)

errs = zero(θ_grid);
for i in eachindex(θ_grid)
    c = [cos((θ_grid[i])*π/180); sin((θ_grid[i])*π/180)]
    p = Φ_sym *c
    c2 = [cos((θ_grid[i])*π/180); sin((θ_grid[i])*π/180)]
    # errs[i] = norm(p./c2)
    errs[i] = (p[1]/(c2[1] + 1e-4) - p[2]/(c2[2] + 1e-4))^2
end

argmin(errs)
θ_rec = θ_grid[argmin(errs)]

=#

a = norm(Φ_obs * [cos((θ1 + β)*π/180); sin((θ1 + β)*π/180)] , 2)
b = norm(Φ_obs * [cos((θ1 + β + 90)*π/180); sin((θ1 + β + 90)*π/180)] , 2)

Δθ = [cos(1* π/180 ) cos(2* π/180 ); sin(1* π/180 ) sin(2* π/180 )]
ps = Φ_sym * Δθ
vec_angle(v::Vector) = angle(v[1] + im * v[2]) * 180/π

vec_angle(ps[:,2]) - vec_angle(Δθ[:,2])
vec_angle(ps[:,1]) - vec_angle(Δθ[:,1])

dangle = vec_angle(ps[:,2]) - vec_angle(ps[:,1])

Φ_rec = R(θ_rec) * Φ_sym * R(-θ_rec) ## Φ_true
Φ_true
Φ_pred = R(β) * R(-θ_rec) * Φ_rec * R(θ_rec) * R(β) ## Φ_obs
Φ_obs

Z_pred = R(θ_rec) * Z_obs * R(-θ_rec) ## Z_true
Z_true




Δθ = [cos(1* π/180 ) cos(2* π/180 ); sin(1* π/180 ) sin(2* π/180 )]
ps = Φ_obs * Δθ
vec_angle(v::Vector) = angle(v[1] + im * v[2]) * 180/π

vec_angle(ps[:,2]) - vec_angle(Δθ[:,2])
vec_angle(ps[:,1]) - vec_angle(Δθ[:,1])

dangle = vec_angle(ps[:,2]) - vec_angle(ps[:,1])

eigen(Φ_sym)





##  ====== 3D =========
x_obs1 = rand(2)
x_obs2 = rand(2) .* 0.1
X_obs = [x_obs2[1] x_obs1[1]; -x_obs1[2] x_obs2[2]] 
y_obs1 = rand(2)
y_obs2 = rand(2) .* 0.1
Y_obs = [y_obs2[1] y_obs1[1]; -y_obs1[2] y_obs2[2]] 

Z_obs = @. X_obs + im * Y_obs

Φ_obs = inv(X_obs) * Y_obs
# Z_21 and -Z_21 yield the same Φ_obs

get_β(M) = 0.5 * atan((M[1,2] - M[2,1])/(M[1,1] + M[2,2])) * 180/π
β = get_β(Φ_obs)
# β = 0 

Φ_sym = R(-β) * Φ_obs * R(-β)
# Φ_sym is the same as Φ_obs
u, k, vt = svd(Φ_sym)
#=
u, vt are identity matrices, diagm(k) is equal to Φ_obs
=#















theta2 = 0.5 * atan((2Φ_sym[1,2])/(Φ_sym[1,1]- Φ_sym[2,2])) * 180/π
Φ_prime = R(theta2) * Φ_sym * R(-theta2)

Z_prime = R(theta2) * R(-β) * Z_obs * R(-β) * R(-theta2)
Z_prime = R(theta2) * Z_obs * R(-theta2)

R(theta2) * R(-β) * Φ_obs * R(-β) * R(-theta2)

X_prime = R(theta2) * R(-β) * X_obs * R(-β) * R(-theta2)
Y_prime = R(theta2) * R(-β) * Y_obs * R(-β) * R(-theta2)

u, k, vt = svd(Φ_sym)

λ, Q = eigen(Φ_sym)
Λ = diagm(λ)

# function get_angle(vs, vc)
#     if vs > 0 && vc > 0
#         asin(vs)
# end

θ = angle(Q[1,1] + im * Q[2,1]) * 180/ π

R(θ) * Φ_sym * R(-θ)
R(90 +- θ) * Φ_sym * R(90 - θ)

Φ_sym






function R(θ_p)
    θ = θ_p *π/180
    reshape([cos(θ), -sin(θ), sin(θ), cos(θ)], 2,2)
end

det_u = det(u)
det_vt = det(vt)

λ, Q = eigen(Φ_sym)
Λ = diagm(λ)
Q * Λ * inv(Q) .- Φ_obs
det_Q = abs(det(Q))
# Φp = det_Q .* Λ
# Rp = Q ./ (det_Q ^ 0.5)



Φ = diagm(randn(2)) .* 10

## == export

extrema([angle(ix + im * iy) for ix in -1:0.02:1, iy in -1:0.02:1])

extrema([atan(iy/ (ix + 1e-2)) for ix in -1:0.02:1, iy in -1:0.02:1])


atan(1) * 180/π
atan(-1) * 180/π

angle(1 + im * 1) * 180/π
angle(-1 + im * 1) * 180/π
angle(1 - im * 1) * 180/π
angle(-1 - im * 1) * 180/π

