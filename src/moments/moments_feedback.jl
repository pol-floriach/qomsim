# Script for time evolution of the moments and variances (SDE system) with feedback
using DifferentialEquations, StaticArrays, Plots, LaTeXStrings, ProgressLogging

# Moment time evolution (deterministic part)
function moments_evolution(u,p,h,t)
    @views γ0, ω, k, η = p
    @views x, p, Vx, Vp, Cxp = u
    dx    = -γ0/2*x + ω*p
    dp    = -γ0/2*p - ω*x
    dVx   = -γ0*Vx + 2*ω*Cxp - 8*k*η*Vx^2  + γ0*(nth+0.5)
    dVp   = -γ0*Vp - 2*ω*Cxp - 8*k*η*Cxp^2 + γ0*(nth+0.5) +2*k
    dCxp  = ω*Vp - ω*Vx - γ0*Cxp - 8*k*η*Vx*Cxp
    SA[dx, dp, dVx, dVp, dCxp]
end

function moments_infogain(u,p,t)
     k, η = p[3:4]
     Vx, Cxp = @view u[3:2:5]
    sqrt(8*k*η)*SA[Vx, Cxp, 0, 0, 0]
end


# Parameters
ω = 2π*1.1
Q = 1e7
γ0 = ω / Q


# Measurement "rate"
g0 = 2π*465
ncav = 1e5
g = g0*sqrt(ncav)
κ = 2π*45e6
Γqba = 4*g^2/κ /1e6
k = Γqba/2
η = 0.9
G = 1

# Constants
ħ = 1.05457182e-34 * 1e9^2 / 1e6
kB = 1.380649e-23 * 1e9^2 / 1e6^2
T = 0#1e-3
const nth = 1/(exp(ħ*ω/(kB*T))-1)
Vx_th = Vp_th = nth + 0.5

# Initial conditions and SDE parameters
p = (γ0, ω, k, η)
u0 = SA[0, 0, Vx_th, Vp_th, 0.0]
tspan = (0.0, 100) # μs

prob= SDEProblem(moments_evolution, moments_infogain, u0, tspan, p)
@time sol = solve(prob,
    # EM(),
    saveat = 0.1/ω,
    #dtmin = 1e-4,
    # dtmax = 0.01,
    dt = 1e-1,
    maxiters = tspan[2]*1e7,
    progress = true
);


# pmeans = plot(sol.t, sol[1,:], label = L"\langle x \rangle", xlabel = "t [μs]", ylabel = "Mean")
# plot!(sol.t, sol[2,:], label = L"\langle p \rangle")

pvars = plot(sol.t, sol[1,:], label = L"\langle x \rangle", xlabel = "t [μs]", ylabel = "Nondimensional")
plot!(sol.t, sol[2,:], label = L"\langle p \rangle")

G = 0
p2 = (p...,G)
prob2 = ODEProblem(covariances_evolution_feedback, u0, tspan, p2)
@time solfeedback = solve(prob2,
    Tsit5(),
    saveat = 0.1/ω,
    #dtmin = 1e-4,
    dtmax = 0.01,
    dt = 1e-1,
    maxiters = tspan[2]*1e3,
    progress = true
);

plot!(solfeedback.t, solfeedback[1,:], label = L"V_x "*", G = $G", xlabel = "t [μs]", ylabel = "Covariance")
plot!(solfeedback.t, solfeedback[2,:], label = L"V_p"*", G = $G")
end
display(pvars)
