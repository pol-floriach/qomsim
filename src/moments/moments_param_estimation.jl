# Script for time evolution of the moments and variances (SDE system)
using DifferentialEquations, StaticArrays, Plots, LaTeXStrings
using DiffEqParamEstim, Optimization, ForwardDiff, OptimizationOptimJL

# Moment time evolution (deterministic part)
function moments_evolution(u,p,t)
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

# Constants
ħ = 1.05457182e-34 * 1e9^2 / 1e6
kB = 1.380649e-23 * 1e9^2 / 1e6^2
T = 300
const nth = 1/(exp(ħ*ω/(kB*T))-1)
Vx_th = Vp_th = nth + 0.5


# Initial conditions and SDE parameters
p = (γ0, ω, k, η)
u0 = SA[0, 0, Vx_th, Vp_th, 0.0]
tspan = (0.0, 100) # μs
prob= SDEProblem(moments_evolution, moments_infogain, u0, tspan, p)


# Generating synthetic Data with the parameters above
using RecursiveArrayTools
t = collect(tspan[1]:0.005/ω:tspan[2])
function generate_data(t)
    sol = solve(prob, SOSRI(),
    saveat = 0.005/ω,
    dt = 1e-2,
    dtmax = 0.1/ω,
    maxiters = tspan[2]*1e7,)
    randomized = VectorOfArray([(sol(t[i]) +0*randn(5)) for i in eachindex(t)])
    data = convert(Array, randomized)
end

aggregate_data = convert(Array, VectorOfArray([generate_data(t) for i in 1:100]))

# PARAMETER ESTIMATION

monte_prob = EnsembleProblem(prob)
# Optimize
obj = build_loss_objective(monte_prob, SOSRI(), 
                           L2Loss(t, aggregate_data, differ_weight = 0.5, data_weight = 0.7),
                           AutoForwardDiff(),
                           verbose = false, trajectories = 100, maxiters = 1e5)
optprob = OptimizationProblem(obj, collect(p), lb = [0.0, 0.0, 0.0, 0.5], ub = [1e-5, 15, 1e-2, 1])

result = solve(optprob, LBFGS())
result.original