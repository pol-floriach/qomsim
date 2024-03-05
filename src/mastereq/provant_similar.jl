using DifferentialEquations, LinearAlgebra

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
H = Diagonal((ω_mech*at*a))
H_nh = Diagonal(H - im*ħ/2* (γ₀*(n+1)*at*a + γ₀*n*a*at))
H_nh_dagger = Diagonal(adjoint(H_nh))


# Jump operators and their rates
J = [a,at]
rates = [γ₀*(n+1),γ₀*n]

ρ0 = zeros(ComplexF64,n,n); ρ0[1] = 1.0


function solve_master!(tspam,ρ0, H, J, rates; kwargs...)
    
    tmp = copy(ρ0)
    dmaster_(t,ρ,dρ) = dmaster_!(dρ,H,J,Jdagger,rates,ρ,tmp)
    integrate_master!(tspan,dmaster_, ρ0; kwargs...)
end


function dmaster_h!(drho, H, J, Jdagger, rates, rho, drho_cache)
    mul!(drho,H,rho,-eltype(rho)(im),zero(eltype(rho)))
    mul!(drho,rho,H,eltype(rho)(im),one(eltype(rho)))
    for i in eachindex(J)
        mul!(drho_cache,J[i],rho,eltype(rho)(rates[i]),zero(eltype(rho)))
        mul!(drho,drho_cache,Jdagger[i],true,true)

        mul!(drho,Jdagger[i],drho_cache,eltype(rho)(-0.5),one(eltype(rho)))

        mul!(drho_cache,rho,Jdagger[i],eltype(rho)(rates[i]),zero(eltype(rho)))
        mul!(drho,drho_cache,J[i],eltype(rho)(-0.5),one(eltype(rho)))
    end
    return drho
end


function integrate_master!(tspan,df, u0, state, dstate; alg = Tsit5(), kwargs...)

    prob = ODEProblem(df,u0, tspan)
    sol = solve(prob, alg; 
        reltol = 1.0e-6, abstol = 1.0e-8,
        save_everystep = false
    )
    return sol
end