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
     k, η = @view p[3:4]
     Vx, Cxp = @view u[3:2:5]
    sqrt(8*k*η)*SA[Vx, Cxp, 0, 0, 0]
end

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
k = 1.5#.5γ0*Vp_th
η = 0.95

# Constants
ħ = 1.05457182e-34 * 1e9^2 / 1e6
kB = 1.380649e-23 * 1e9^2 / 1e6^2
T = 0#1e-3
const nth = 1/(exp(ħ*ω/(kB*T))-1)
Vx_th = Vp_th = nth + 0.5

# Initial conditions and SDE parameters
p = (γ0, ω, k, η)
u0 = SA[Vx_th, Vp_th, 0.0]
tspan = (0.0, 10) # μs

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

#plot(sol.t, sol[3,:])
# Data visualization
# xlabels = [L"\langle x \rangle", L"\langle p \rangle", L"V_x", L"V_p", L"C_{xp}"];
# xlabels = [L"V_x", L"V_p"];
# myplot2 = [plot(sol.t, (sol[i,:]), ylabel = xlabels[i]) for i in 1:2];
# plot(myplot2..., 
#     layout = (2,1),
#     dpi = 1400, 
#     legend = false, 
#     #yscale = :log10,
#     #plot_title = "Moments time evolution, \n with thermal initial state",
#     #size = (600,700),
# )
# pmeans = plot(sol.t, sol[1,:], label = L"\langle x \rangle", xlabel = "t [μs]", ylabel = "Mean")
# plot!(sol.t, sol[2,:], label = L"\langle p \rangle")

pvars = plot(sol.t, sol[1,:], label = L"V_x", xlabel = "t [μs]", ylabel = "Covariance")
plot!(sol.t, sol[2,:], label = L"V_p")
# plot(pmeans, pvars, layout = (2,1))
