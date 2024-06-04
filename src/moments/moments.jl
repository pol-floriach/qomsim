# Script for time evolution of the moments and variances (SDE system)
using DifferentialEquations, StaticArrays, Plots, LaTeXStrings, ProgressLogging

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
η = 0.0#0.9
G = 1

# Constants
ħ = 1.05457182e-34 * 1e9^2 / 1e6
kB = 1.380649e-23 * 1e9^2 / 1e6^2
T = 300
const nth = 1/(exp(ħ*ω/(kB*T))-1)
Vx_th = Vp_th = nth + 0.5

# Initial conditions and SDE parameters
p = (γ0, ω, k, η)
u0 = SA[1.0, 0.0, Vx_th, Vp_th, 0.0]
tspan = (0.0, 100) # μs

# prob= SDEProblem(moments_evolution, moments_infogain, u0, tspan, p)
prob = ODEProblem(moments_evolution, u0, tspan, p)
@time sol = solve(prob,
    # EM(),
    saveat = 0.1/ω,
    #dtmin = 1e-4,
    # dtmax = 0.01,
    dt = 1e-1,
    maxiters = tspan[2]*1e7,
    progress = true,
    save_idxs = 1
);


# pmeans = plot(sol.t, sol[1,:], label = L"\langle x \rangle", xlabel = "t [μs]", ylabel = "Mean")
# plot!(sol.t, sol[2,:], label = L"\langle p \rangle")

pvars = plot(sol.t, sol[1,:], label = L"\langle x \rangle", xlabel = "t [μs]", ylabel = "Nondimensional")
plot!(sol.t, sol[2,:], label = L"\langle p \rangle")



# Parameter estimation for ODE
using DiffEqParamEstim, Optimization, ForwardDiff, OptimizationOptimJL
t = tspan[1]:0.1/ω:tspan[2]
data = sol.u + 0.0001*randn(size(sol))

cost_function = build_loss_objective(prob, Tsit5(), L2Loss(t, data),
                                    AutoForwardDiff(),
                                     maxiters = 10000, verbose = false, save_idxs = 1
)


optprob = OptimizationProblem(cost_function, collect(p), lb = [0.0, 0.0, -0.5, -0.5], ub = [1e-6, 15.0, 0.1, 1.0])

result = solve(optprob, BFGS())
result.original
result.u