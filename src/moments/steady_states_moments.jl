# Script that computes steady state values of covariances of the quadratures
using DifferentialEquations, StaticArrays, Plots, LaTeXStrings
function covariances_ss(u,p)
    @views γ0, ω, k, η, nth = p
    @views Vx, Vp, Cxp = u
    dVx   = -γ0*Vx + 2*ω*Cxp - 8*k*η*Vx^2  + γ0*(nth+0.5)   
    dVp   = -γ0*Vp - 2*ω*Cxp - 8*k*η*Cxp^2 + γ0*(nth+0.5) + 2*k  
    dCxp  = ω*Vp - ω*Vx - γ0*Cxp - 8*k*η*Vx*Cxp
    SA[dVx, dVp, dCxp]
end

# Parameters 
m = 1e-12   # [m] = kg
ω = 2π*1.1  # [ω] = MHz
Q = 1e9
γ0 = ω / Q  # [γ0] = MHz



# Constants
const ħ = 1.05457182e-34 * 1e9^2 / 1e6    # [ħ] = kg⋅nm^2⋅μs^-1
const kB = 1.380649e-23 * 1e9^2 / 1e6^2   # [kB] = kg⋅nm^2⋅μs^-2

η = 0.8
T = 1e2


nth = 1/(exp(ħ*ω/(kB*T))-1)
xzp = sqrt(ħ/(2*m*ω))/1e9               # [x_zp] = m

# Steady state / initial condition for continuation starting at k̄ = 0
Vx_th = Vp_th = nth + 0.5                # Unitless
# Initializing vectors
kvec = collect(0:1e-1:10)
steady_values = zeros(length(kvec), 3)
snrvec = zeros(length(kvec))
u0_newt = SA[Vx_th, Vp_th, 0.0]

for ii in eachindex(kvec)
    # Update parameter tuple
    p = (γ0, ω, kvec[ii], η, nth)
    # Find steady state
    prob_ss = NonlinearProblem(covariances_ss, u0_newt, p)
    ss = solve(prob_ss, NewtonRaphson(), reltol = 1e-10, abstol = 1e-10)
    steady_values[ii,:] = ss.u
    # Update initial guess for next Newton iteration
    u0_newt = ss.u
    # Find SNR
    # snrvec[ii] = (abs(ss[1]*ħ/(m*ω) - ss[2]*m*ħ*ω))*2*η*kvec[ii]*m*ω/ħ
    snrvec[ii] = abs(ss[1] - ss[2])*2*η*kvec[ii]
end

# Plot of Variances
vars = hline([0.5], label = L"\frac{1}{2}")
plot!(kvec, steady_values[:,1], 
    label = L"V_x",
    xlabel = L"\tilde k", 
    ylabel = "Variances (nondimensional)",
    lt = :scatter,
    ms = 1.5,
    msw = 0,
    ma = 0.5,
)
plot!(kvec, steady_values[:,2],
    label = L"V_p",
    lt = :scatter,
    ms = 1.5,
    msw = 0,
    ma = 0.5,
    ylims = (0,2)
)
#vline!([4*xzp^2/γ0])

# Product of variances (nondimensional) vs heisenberg unc. lower bound
product = plot(kvec[2:end], steady_values[2:end,1].*steady_values[2:end,2], label = L"V_x\cdot V_p",ylabel = "Product of variances",ylims = (0,1))
hline!([0.25], label = L"\frac{1}{4}")#, ylims = (0,1))

# "Transduction" noise vs zero point motion uncertainty? TODO: not only ground state
SN = [1/(η*k*1e6)*xzp^2 for k in kvec]
Sxx = 8*(nth+0.5)*xzp^2/(γ0*1e6)
Sxx_gs = 4*xzp^2/(γ0*1e6)

myplot = plot(kvec[2:end], SN[2:end], label = L"S_N = \frac{1}{2\eta k}", yscale = :log10, xlabel = L"\tilde k",ylabel = "Noise [m²/Hz]")#, ylims = (1e-30, 1e-21))
hline!([Sxx], label = L"S_{xx} = \frac{8x_{zp}^2}{\gamma_0}(\bar n + \frac{1}{2})")
hline!([Sxx_gs], label = L"S_{xx}^{gs}")

# Signal to noise ratio vs k̃
snrplot = plot(kvec[2:end], snrvec[2:end],
    # yscale = :log10,
    legend = false,
    ylabel = "SNR",
    title = "SNR = "*L"\frac{|V_x - V_p|}{S_N}"*" for η = $(η)", ylims = (-1,25)
)
hline!([1])

total = plot(vars, myplot, product, snrplot, layout = (2,2), size = (900,900), xlabel = L"\tilde k", plot_title = "η = $(η), T = $(T), Q = $(Q)")
display(total)


# Finding out what the measurement rate is!

g0 = 2π*465
ncav = 1e5
g = g0*sqrt(ncav)
κ = 2π*45e6
Γqba = 4*g^2/κ /1e6# [Γqba] = MHz

Sxximp = 1/(2*Γqba)* xzp^2

ωl = 2π* 3e8 / 1.55e-6
ħ2 = 1.05e-34

# Checking out the laser power needed to have a specifc measurement rate
P = ncav * ħ2 * ωl * (κ/2)^2 / κ


kvec2 = 0:1e-8:1e-6
SN = [1/(η*k*1e6)*xzp^2 for k in kvec2]
Sxx = 8*(nth+0.5)*xzp^2/(γ0*1e6)
Sxx_gs = 4*xzp^2/(γ0*1e6)

myplot2 = plot(kvec2, SN, label = L"S_N = \frac{1}{2\eta k}", yscale = :log10, xlabel = L"\tilde k",ylabel = "Noise [m²/Hz]")#, ylims = (1e-30, 1e-21))
hline!([Sxx], label = L"S_{xx} = \frac{8x_{zp}^2}{\gamma_0}(\bar n + \frac{1}{2})")
hline!([Sxx_gs], label = L"S_{xx}^{gs}")