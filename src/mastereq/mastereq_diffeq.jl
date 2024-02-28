using DifferentialEquations, Plots, LinearAlgebra

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

function master_eq!(dρ,ρ,param,t)
    a, at, H, n,  ħ, γ₀ = param
    @views
     dρ .= -im/ħ * (H*ρ - ρ*H) + γ₀*(n+1)*(a*ρ*at - at*a*ρ/2 - ρ/2*at*a) + γ₀*n*(at*ρ*a - a*at*ρ/2 - ρ/2*a*at)
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
a = anihilation(n);
at = creation(n);
H = ω_mech*at*a;

p = (a, at, H, n, ħ, γ₀)

# Initial condition (ground state)
ρ0 = zeros(ComplexF64,n,n); ρ0[1,1] = 1;

# Simulation parameters
tspan = (0.0,1e6)

prob = ODEProblem(master_eq!, ρ0, tspan, p);

@time sol = solve(prob,Euler(),saveat = 1000, dt = 1e-1, dtmin = 1e-2, dtmax = 1e-1, progress = true);

expect = expectation(at*a,sol)

# heatmap(abs2.(sol[end]))

plot(0:1000:1e6,expect, xlabel = L"t [\mu s]", ylabel = L"\langle n_m \rangle", legend = false, dpi = 1400)