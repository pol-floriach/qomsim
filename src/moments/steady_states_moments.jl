using DifferentialEquations, StaticArrays, Plots, LaTeXStrings
function covariances_ss(u,p)
    @views γ0, ω, k, η = p
    @views Vx, Vp, Cxp = u
    dVx   = -γ0*Vx + 2*ω*Cxp - 8*k*η*Vx^2  + γ0*(nth+0.5)   
    dVp   = -γ0*Vp - 2*ω*Cxp - 8*k*η*Cxp^2 + γ0*(nth+0.5) + 2*k  
    dCxp  = ω*Vp - ω*Vx - γ0*Cxp - 8*k*η*Vx*Cxp
    SA[dVx, dVp, dCxp]
end

function snrmeasure(sol,η,k, xzp)
    SNR = (abs(sol[1] - sol[2]))*(η*k)/xzp^2
end

# Parameters
m = 1e-12
ω = 2π*1.1
Q = 1e7
γ0 = ω / Q
η = 0.9


# Constants
ħ = 1.05457182e-34 * 1e9^2 / 1e6
kB = 1.380649e-23 * 1e9^2 / 1e6^2
T = 0#1e-3
const nth = 1/(exp(ħ*ω/(kB*T))-1)
Vx_th = Vp_th = nth + 0.5

xzp = sqrt(ħ/(2*m*ω))/1e9 # en m


kvec = collect(0:1e-1:10)
steady_values = zeros(length(kvec), 3)
snrvec = zeros(length(kvec))
ii = 0
u0_newt = SA[Vx_th, Vp_th, 0.0]

for ii in eachindex(kvec)
    # Update parameter tuple
    p = (γ0, ω, kvec[ii], η)
    # Find steady state
    prob_ss = NonlinearProblem(covariances_ss, u0_newt, p)
    ss = solve(prob_ss, NewtonRaphson(), reltol = 1e-10, abstol = 1e-10)
    steady_values[ii,:] = ss.u
    # Update initial guess for next Newton iteration
    u0_newt = ss.u
    # Find SNR
    snrvec[k] = snrmeasure(ss,p, k)
end

vars = hline([0.5], label = L"\frac{1}{2}")
plot!(kvec, steady_values[:,1], 
    label = L"V_x",
    xlabel = L"\tilde k", 
    ylabel = "Variances (nondimensional)",
    lt = :scatter,
    ms = 1.5,
    msw = 0,
    ma = 0.5,
    title = "η = $(η)"
)
plot!(kvec, steady_values[:,2],
    label = L"V_p",
    lt = :scatter,
    ms = 1.5,
    msw = 0,
    ma = 0.5
)
vline!([4*xzp^2/γ0])

product = plot(kvec, steady_values[:,1].*steady_values[:,2], label = L"V_x\cdot V_p",ylabel = "Product of variances",)
hline!([0.25], label = L"\frac{1}{4}", ylims = (0,1))

# Soroll?
SN = [1/(η*k)*xzp^2 for k in kvec]
Sxx = xzp^2/(γ0*1e6)

myplot = plot(kvec[2:end], SN[2:end], label = L"S_N = \frac{1}{2\eta k}", yscale = :log10, xlabel = L"\tilde k",ylabel = "Noise [m²/Hz]")
hline!([8*nth*Sxx], label = L"S_{xx} = \frac{x_{zp}^2}{\gamma_0}")

plot(vars, myplot, product,layout = (3,1), size = (700,700), xlabel = L"\tilde k")
