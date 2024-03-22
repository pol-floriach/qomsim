# Mechanical oscillator in heat bath

using QuantumOptics, Plots, LaTeXStrings, ProgressLogging, DifferentialEquations

# Parameters
ω_mech = 2π*1.1
Q = 1e7
γ0 = ω_mech / Q

# Constants
ħ = 1.05457182e-34 * 1e9^2 / 1e6
kB = 1.380649e-23 * 1e9^2 / 1e6^2
T = 1e-3
n = 1/(exp(ħ*ω_mech/(kB*T))-1)

# Hilbert space
b_mech = FockBasis(100)

# Operators
b = destroy(b_mech)
bt = create(b_mech)
bt = dagger(b)
N = number(b_mech)

# Mechanical oscillator Hamiltonian
H_mech = ω_mech*bt*b
#H_mech2 = ω_mech*N


# Constants for k
m = 1e-12
xzpf = sqrt(ħ/(2*m*ω_mech))
nth = 1/(exp(ħ*ω_mech/(kB*T))-1)

k = γ0 * nth #/ xzpf^2

x = b+bt


J = [b,bt]
rates = [γ0*(n+1),γ0*n]


# J2 = [b,bt,x]
# rates2 = [γ0*(n+1),γ0*n,k]

J2 = [x]
rates2 = [k]

# Initial conditions
Ψ0 = fockstate(b_mech,0)
ρ0 = dm(Ψ0)


using LinearAlgebra
BLAS.set_num_threads(1)
tspan = collect(0:1e3:1e6)
dt = 0.1
@time tout, ρt = timeevolution.master(tspan,ρ0, H_mech, J; dt = dt, rates = rates, progress= true);
@time tout2, ρt2 = timeevolution.master(tspan,ρ0, H_mech, J2; dt = dt, rates = rates2, progress= true);


plot(tspan,real(expect(bt*b,ρt)),
    title = "Mechanical oscillator",
    xlabel = L"t [\mu s]",
    ylabel = L"\langle n_m \rangle",
    label = "w/o BA",
    dpi = 1400,
)
# plot!(tout2,real(expect(bt*b,ρt2)),
#     xlabel = L"t [\mu s]",
#     ylabel = L"\langle n_m \rangle",
#     label = "w/ BA",
#     #lc = 2,
#     dpi = 1400,
#  #   xlims = (0,*2π/ω_mech)
# )

# Solucio analitica (potser)



Nop = Diagonal(Matrix(N.data))
G = 2*nth*sinh.(γ0/2*tspan) ./ (cosh.(γ0/2*tspan) + (1+2*nth)*sinh.(γ0/2*tspan))
G = G[2:end]

ρex = [(1 .- G[i])*exp(log.(G[i])*Nop) for i in eachindex(G)]


expected_n = [real(tr(ρex[i]*Nop)) for i in eachindex(G)]


plot!(tspan[2:10:end], expected_n[1:10:end],
    label = "Analytical",
    dpi = 1400,
)




# # Stochastic master eqn.
# J = [b,bt,1/sqrt(2)*(b+bt)]
# rates = [γ0*(n+1),γ0*n,γ0]
# C = [sqrt(γ0)/sqrt(2)*(b+bt)]
# C = [0*(b+bt)]
# dt = 1e-3
# tspan = collect(0:10:1e4)
# @time tout, ρt_stoch = stochastic.master(tspan, ρ0, H_mech, J, C; dt = dt, rates = rates, alg = EM(), progress = true);

# plot(tout,real(expect(bt*b,ρt_stoch)),
#     dpi = 1400,
#     lt = :scatter,
#     msw = 0,
#     ms = 1,
#     ylims = (0,0.2)
# )

# heatmap(real(ρt_stoch[end].data), xlims = (0,10),ylims = (0,10)),