using DifferentialEquations, Plots, LinearAlgebra, ProgressLogging, LaTeXStrings, SparseArrays

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

# With backaction
function master_eq_mul!(dρ, ρ, param, t)
    a, at, x, H, γ₀n1, γ₀n, k, γiJρ = param
    mul!(dρ,H,ρ,-eltype(ρ)(im),zero(eltype(ρ)))   # im/ħ ρ⋅H 
    mul!(dρ,ρ,H,eltype(ρ)(im),one(eltype(ρ)))     # -im/ħ H⋅ρ

    # Terme a
    mul!(γiJρ,a,ρ,eltype(ρ)(γ₀n1),zero(eltype(ρ)))   # γᵢaρ
    mul!(dρ,γiJρ,at,true,true)                       # dρ += γᵢa⋅ρ⋅at

    mul!(dρ,at,γiJρ,eltype(ρ)(-0.5),one(eltype(ρ))) # dρ -= 1/2 * γᵢ at⋅a⋅ρ

    mul!(γiJρ,ρ,at,eltype(ρ)(γ₀n1),zero(eltype(ρ))) # γᵢ ρ⋅at
    mul!(dρ,γiJρ,a,eltype(ρ)(-0.5),one(eltype(ρ)))  # dρ -= 1/2 * γᵢ ρ⋅at⋅a

    # Terme at
    mul!(γiJρ,at,ρ,eltype(ρ)(γ₀n),zero(eltype(ρ)))  # γᵢatρ
    mul!(dρ,γiJρ,a,true,true)                       # dρ += γᵢat⋅ρ⋅a

    mul!(dρ,a,γiJρ,eltype(ρ)(-0.5),one(eltype(ρ)))  # dρ -= 1/2 * γᵢ a⋅at⋅ρ

    mul!(γiJρ,ρ,a,eltype(ρ)(γ₀n),zero(eltype(ρ)))   # γᵢ ρ⋅at
    mul!(dρ,γiJρ,at,eltype(ρ)(-0.5),one(eltype(ρ))) # dρ -= 1/2 * γᵢ ρ⋅a⋅at 

    # Terme x
    mul!(γiJρ,x,ρ,eltype(ρ)(k),zero(eltype(ρ)))     # k⋅xρ
    mul!(dρ,γiJρ,x,true,true)                       # dρ += k⋅x⋅ρ⋅xt

    mul!(dρ,x,γiJρ,eltype(ρ)(-0.5),one(eltype(ρ)))  # dρ -= 1/2 * k⋅xt⋅x⋅ρ

    mul!(γiJρ,ρ,x,eltype(ρ)(γ₀n),zero(eltype(ρ)))   # k⋅ρ⋅at
    mul!(dρ,γiJρ,x,eltype(ρ)(-0.5),one(eltype(ρ))) # dρ -= 1/2 * k⋅ρ⋅a⋅at 
end


function master_eq_nh!(dρ, ρ, param,t)
    a, at, H_nh, ħ = param
    dρ .= -im/ħ * (H_nh*ρ - ρ*H_nh) + a*ρ*at + at*ρ*a
end

function master_eq_nh_mul!(dρ, ρ, param, t)
    a, at, H_nh, ħ, γ₀n1, γ₀n, γiJρ = param
    mul!(dρ, H_nh, ρ, -eltype(ρ)(im)/ħ, zero(eltype(ρ)))
    mul!(dρ, ρ, H_nh, eltype(ρ)(im)/ħ, one(eltype(ρ)))

    # γ₀(n+1) term
    mul!(γiJρ,a,ρ,eltype(ρ)(γ₀n1),zero(eltype(ρ)))   # γᵢaρ
    mul!(dρ,γiJρ,at,true,true)                       # dρ += γᵢa⋅ρ⋅at

    # γ₀⋅n term
    mul!(γiJρ,at,ρ,eltype(ρ)(γ₀n),zero(eltype(ρ)))  # γᵢatρ
    mul!(dρ,γiJρ,a,true,true)                       # dρ += γᵢat⋅ρ⋅a
end

function expectation(operator, ρ)
    N = size(ρ,3)
    mean_operator = [real(tr(ρ[:,:,i]*operator)) for i in 1:N]
end

# Parameters
ω_mech = 2π*1.1
Q = 1e7
γ₀ = ω_mech / Q

# Constants
ħ = 1.05457182e-34 * 1e9^2 / 1e6
kB = 1.380649e-23 * 1e9^2 / 1e6^2
T = 1e-3
n = 1/(exp(ħ*ω_mech/(kB*T))-1)

N = 50
# Operators
a = Bidiagonal(anihilation(N+1), :U)
at = Bidiagonal(creation(N+1),:L)
H = Diagonal(ω_mech*at*a)
# H_nh = Diagonal(H - im*ħ/2* (γ₀*(n+1)*at*a + γ₀*n*a*at))
# H_nh_dagger = Diagonal(adjoint(H_nh))



x = a+at
k = γ₀*n

# Initial condition (ground state)

# Initial conditions for coherent state
# ψ0 = zeros(ComplexF64,N+1)
# α = 3
# for i in 0:20
#     ψ0[i+1] = exp(-1/2*abs2(α))*α^i/sqrt(factorial(i))
# end
# ρ0 = ψ0 * adjoint(ψ0)

ρ0 = zeros(ComplexF64,N+1,N+1); ρ0[1,1] = 1;

# p = (a, at, H, n, ħ, γ₀)
# p2 = (a,at, H_nh, ħ)

# amb backaction
param =  (a, at, x, H, γ₀*(n+1), γ₀*n, k, similar(ρ0))

# param2 = (a, at, H_nh, H_nh_dagger, ħ, γ₀*(n+1), γ₀*n, similar(ρ0))

# Simulation parameters
tspan = (0.0,1e4)

prob = ODEProblem(master_eq_mul!, ρ0, tspan, param)
# prob = ODEProblem(master_eq!, ρ0, tspan, p)
@time sol = solve(prob,
    RK4(),
    saveat = 1e1, 
    # dtmax = 1e-2, 
    dtmin = 1e-3,
    dt = 1e-1, 
    progress = true;
    reltol = 1e-5, abstol = 1e-7,
    save_everystep=false,
    maxiters = tspan[2]*1e3
);

expect = expectation(at*a,sol)

plot!(sol.t,expect,
        xlabel = L"t [\mu s]", 
        ylabel = L"\langle n_m \rangle", 
        legend = false, 
        dpi = 1400,
)




# prob2 = ODEProblem(master_eq_nh!, ρ0, tspan, p2)
# @time sol2 = solve(prob, Tsit5(), saveat = 1, dt = 1e-1, dtmin = 1e-2, dtmax = 1e-1, progress = true, maxiters);
# expect2 = expectation(at*a,sol2)

# plot(0:1000:1e5, expect2)


# prob3 = ODEProblem(master_eq_nh_mul!, ρ0, tspan, param2)
# @time sol = solve(prob,
#                 Tsit5(),
#                 saveat = 10, 
#                 dtmax = 1e-2, 
#                 dtmin = 1e-3,
#                 dt = 1e-1, 
#                 progress = true;
# );

# expect3 = expectation(at*a,sol)
# plot(0:10:1000-1, expect3)
