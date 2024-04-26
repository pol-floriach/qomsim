using DifferentialEquations, StaticArrays, Plots

# Moment time evolution (deterministic part)
function moments_evolution(u,p,t)
    @views γ0, ω, k, η = p
    @views x, p, Vx, Vp, Cxp = u
    dx    = -γ0/2*x + ω*p
    dp    = -γ0/2*p - ω*x
    dVx   = -γ0*Vx + 2*ω*Cxp - 8*k*η*Vx^2  + γ0*(nth+0.5)
    dVp   = -γ0*Vp - 2*ω*Cxp - 8*k*η*Cxp^2 + γ0*(nth+0.5) +2*k
    dCxp  = ω*Vp - ω*Vx - γ0*Cxp - 8*k*η*Vx*Cxp
end

function moments_infogain(u,p,t)
    sqrt(8*k*η)*SA[Vx, Cxp, 0, 0, 0]
end

# Parameters
ω = 2π*1.1
Q = 1e7
γ0 = ω / Q
k = 0.0
η = 0.0

# Constants
ħ = 1.05457182e-34 * 1e9^2 / 1e6
kB = 1.380649e-23 * 1e9^2 / 1e6^2
T = 1e-3
const nth = 1/(exp(ħ*ω/(kB*T))-1)


# Hilbert space
bb = exp(-ħ*ω/(kB*T))
Vx_th = Vp_th = 1/(1-bb) * (0.5*1/(1-bb) + bb/(1-bb)^2)


# Initial conditions and SDE parameters
p = (γ0, ω, k, η)
u0 = SA[0.0, 0.0, Vx_th, Vp_th, 0.0]
tspan = (0.0, 1e2) # μs

prob= ODEProblem(moments_evolution, u0, tspan, p)
@time sol = solve(prob,
    Tsit5(),
    saveat = 1e-1,
    dtmin = 1e-3,
    dt = 1e-1,
    maxiters = tspan[2]*1e3
)

# Data visualization
xlabels = ["x", "p", "V_x", "V_p", "C_xp"]
myplot = [plot(sol.t, sol[i,:], ylabel = xlabels[i]) for i in 1:5];
plot(myplot..., layout = (5,1), dpi = 1400)
plot(px, pp, pVx, pVp, pCxp, layout = (5,1), dpi = 1400)