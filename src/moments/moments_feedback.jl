# Script for time evolution of the moments and variances (SDE system) with feedback
using StochasticDelayDiffEq, StaticArrays, Plots, LaTeXStrings, ProgressLogging

# Moments time evolution (deterministic part) with feedback term in Hamiltonian: H_fb(t) = -G⋅⟨x(t-T/4)⟩x(t)
function moments_evolution(u,h,p,t)
    @views γ0, ω, k, η, G, t0 = p
    @views x, p, Vx, Vp, Cxp = u
    dx    = -γ0/2*x + ω*p
    dp    = -γ0/2*p - ω*x + G*h(p,t-t0)[1]
    dVx   = -γ0*Vx + 2*ω*Cxp - 8*k*η*Vx^2  + γ0*(nth+0.5)
    dVp   = -γ0*Vp - 2*ω*Cxp - 8*k*η*Cxp^2 + γ0*(nth+0.5) +2*k
    dCxp  = ω*Vp - ω*Vx - γ0*Cxp - 8*k*η*Vx*Cxp
    SA[dx, dp, dVx, dVp, dCxp]
end

# Stochastic part of the equations (due to information gain)
function moments_infogain(u,h,p,t)
     k, η = p[3:4]
     Vx, Cxp = @view u[3:2:5]
    sqrt(8*k*η)*SA[Vx, Cxp, 0, 0, 0]
end


# Parameters
ω = 2π*1.1 # [ω] = MHz
Q = 1e7
γ0 = ω / Q
G = 0
t0 = 2π/ω * 0.25

# Measurement "rate" ([MHz])
g0 = 2π*465
ncav = 1e5
g = g0*sqrt(ncav)
κ = 2π*45e6
Γqba = 4*g^2/κ /1e6
k = Γqba/2
η = 0.9
G = 1

# Constants
const ħ = 1.05457182e-34 * 1e9^2 / 1e6
const kB = 1.380649e-23 * 1e9^2 / 1e6^2
T = 0
const nth = 1/(exp(ħ*ω/(kB*T))-1)
Vx_th = Vp_th = nth + 0.5

# Initialize feedback at zero
h(p,t) = @SVector zeros(5)
lags = [t0]

# Initial conditions and SDE parameters
p = (γ0, ω, k, η, G, t0)
u0 = [0, 0, Vx_th, Vp_th, 0.0]
tspan = (0.0, 1000) # μs

# Simulation
prob= SDDEProblem(moments_evolution, moments_infogain, u0, h, tspan, p)#; constant_lags = lags)
@time sol = solve(prob,
    SOSRI(),
    saveat = 0.1/ω,
    # dtmin = 1e-4,
    # dtmax = 0.01,
    dt = 1e-1,
    maxiters = tspan[2]*1e7,
    progress = true
);

pvars = plot(sol.t, sol[1,:], label = L"\langle x \rangle", xlabel = "t [μs]", ylabel = "Nondimensional")
plot!(sol.t, sol[2,:], label = L"\langle p \rangle", title = "Feedback = -G· "*L"\langle x(t-T/4) \rangle")



# # Provar de fer ensemble problem
# ensembleprob = EnsembleProblem(prob)
# sol = solve(ensembleprob, SOSRI(), saveat = 0.1/ω, dt = 1e-1, maxiters = tspan[2]*1e7, progress = true, trajectories = 10)
# summ = EnsembleSummary(sol)#, 0:0.1/ω:tspan[2])
# plot(summ)