using DifferentialEquations, Plots, LinearAlgebra, ProgressLogging, LaTeXStrings

function anihilation(n)
    a = zeros(n,n)
    for i in 1:n-1
        a[i,i+1] = sqrt(i)
    end
    return a
end

function creation(n)
    at = zeros(n,n)
    for i in 1:n-1
        at[i+1,i] = sqrt(i)
    end
    return at
end

function master_eq!(dρ,ρ,p,t)
    a, at, H, n,  ħ, γ₀ = p
    dρ .= -im/ħ * (H*ρ - ρ*H) + γ₀*(n+1)*(a*ρ*at - at*a*ρ/2 - ρ/2*at*a) + γ₀*n*(at*ρ*a - a*at*ρ/2 - ρ/2*a*at)
end

function master_eq_mul!(dρ, ρ, param, t)
    a, at, H, γ₀n1, γ₀n, γiJρ = param
    mul!(dρ,H,ρ,-eltype(ρ)(im)/ħ,zero(eltype(ρ)))   # im/ħ ρ⋅H 
    mul!(dρ,ρ,H,eltype(ρ)(im)/ħ,one(eltype(ρ)))     # -im/ħ H⋅ρ

    # Terme γ⋅(n+1)
    mul!(γiJρ,a,ρ,eltype(ρ)(γ₀n1),zero(eltype(ρ)))   # γᵢaρ
    mul!(dρ,γiJρ,at,true,true)                       # dρ += γᵢa⋅ρ⋅at

    mul!(dρ,at,γiJρ,eltype(ρ)(-0.5),one(eltype(ρ))) # dρ -= 1/2 * γᵢ at⋅a⋅ρ

    mul!(γiJρ,ρ,at,eltype(ρ)(γ₀n1),zero(eltype(ρ))) # γᵢ ρ⋅at
    mul!(dρ,γiJρ,a,eltype(ρ)(-0.5),one(eltype(ρ)))  # dρ -= 1/2 * γᵢ ρ⋅at⋅a

    # Terme γ⋅n
    mul!(γiJρ,at,ρ,eltype(ρ)(γ₀n),zero(eltype(ρ)))  # γᵢatρ
    mul!(dρ,γiJρ,a,true,true)                       # dρ += γᵢat⋅ρ⋅a

    mul!(dρ,a,γiJρ,eltype(ρ)(-0.5),one(eltype(ρ)))  # dρ -= 1/2 * γᵢ a⋅at⋅ρ

    mul!(γiJρ,ρ,a,eltype(ρ)(γ₀n),zero(eltype(ρ)))   # γᵢ ρ⋅at
    mul!(dρ,γiJρ,at,eltype(ρ)(-0.5),one(eltype(ρ))) # dρ -= 1/2 * γᵢ ρ⋅a⋅at 
end



function master_eq_nh!(dρ, ρ, param,t)
    a, at, H_nh, ħ = param
    dρ .= -im/ħ * (H_nh*ρ - ρ*H_nh) + a*ρ*at + at*ρ*a
end

function expectation(operator, ρ)
    N = size(sol,3)
    mean_operator = zeros(N)
    for i in 1:N
        mean_operator[i] = real(tr(ρ[:,:,i]*operator))
    end
    return mean_operator
end

# Parameters
ω_mech = 2π*1.1
Q = 1e7
γ₀ = ω_mech / Q

# Constants
ħ = 1.05457182e-34 * 1e9^2 / 1e6
kB = 1.380649e-23 * 1e9^2 / 1e6^2
T = 1e-3
n = Int(round(1/(exp(ħ*ω_mech/(kB*T))-1)))


# Operators
a = anihilation(n)
at = creation(n)
H = Diagonal(ω_mech*at*a)
H_nh = Diagonal(H - im*ħ/2* (a*at + at*a))



# Initial condition (ground state)
ρ0 = zeros(ComplexF64,n,n); ρ0[1,1] = 1;

p = (a, at, H, n, ħ, γ₀)
p2 = (a,at, H_nh, ħ)

param = (a, at, H, γ₀*(n+1), γ₀*n, similar(ρ0))

# Simulation parameters
tspan = (0.0,1e3)

prob = ODEProblem(master_eq_mul!, ρ0, tspan, param)
# prob = ODEProblem(master_eq!, ρ0, tspan, p)
@time sol = solve(prob,
                DP5(),
                saveat = 1, 
                dtmax = 1e-1, 
                dtmin = 1e-3,
                dt = 1e-1, 
                progress = true;
);

expect = expectation(at*a,sol)

# heatmap(abs2.(sol[end]))

# plot(0:1000:1e6,expect, xlabel = L"t [\mu s]", ylabel = L"\langle n_m \rangle", legend = false, dpi = 1400)


prob2 = ODEProblem(master_eq_nh!, ρ0, tspan, p2)
@time sol2 = solve(prob, Euler(), saveat = 1000, dt = 1e-1, dtmin = 1e-2, dtmax = 1e-1, progress = true);
expect2 = expectation(at*a,sol2)



function hola(params,p2)
    ρ = rand(ComplexF64, 18,18)
    dρ1 = similar(ρ)
    dρ2 = similar(ρ)
    dρ3 = similar(ρ)

    dmaster_h!(dρ1, H, a, at, ρ, params)
    master_eq_mul!(dρ2, ρ, params, 0.0)
    master_eq!(dρ3,ρ,p2,0.0)
    return dρ1, dρ2, dρ3
end