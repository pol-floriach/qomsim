# Script for time evolution of variances only
using DifferentialEquations, StaticArrays, Plots, LaTeXStrings, ProgressLogging

function covariances_evolution(u,p,t)
    γ0, ω, k, η = p
    Vx, Vp, Cxp = u
    dVx   = -γ0*Vx + 2*ω*Cxp + γ0*(nth+0.5) - 8*k*η*Vx^2
    dVp   = -γ0*Vp - 2*ω*Cxp + γ0*(nth+0.5) + 2*k - 8*k*η*Cxp^2
    dCxp  = ω*Vp - ω*Vx - γ0*Cxp - 8*k*η*Vx*Cxp
    SA[dVx, dVp, dCxp]
end


# Parameters
ω = 2π*1.1
Q = 1e7
γ0 = ω / Q

k = 0.01207
η = 0.9


# Constants
ħ = 1.05457182e-34 * 1e9^2 / 1e6
kB = 1.380649e-23 * 1e9^2 / 1e6^2
T = 0#1e-3
const nth = 1/(exp(ħ*ω/(kB*T))-1)
Vx_th = Vp_th = nth + 0.5

# Initial conditions and SDE parameters
p = (γ0, ω, k, η)
u0 = SA[Vx_th, Vp_th, 0.0]
tspan = (0.0, 500) # μs

prob= ODEProblem(covariances_evolution, u0, tspan, p)
@time sol = solve(prob,
    Tsit5(),
    saveat = 0.1/ω,
    #dtmin = 1e-4,
    dtmax = 0.01,
    dt = 1e-1,
    maxiters = tspan[2]*1e3,
    progress = true
);

# pmeans = plot(sol.t, sol[1,:], label = L"\langle x \rangle", xlabel = "t [μs]", ylabel = "Mean")
# plot!(sol.t, sol[2,:], label = L"\langle p \rangle")

pvars = plot(sol.t, sol[1,:], label = L"V_x", xlabel = "t [μs]", ylabel = "Covariance")
plot!(sol.t, sol[2,:], label = L"V_p", legend = :outertopright, size = (800,400))
