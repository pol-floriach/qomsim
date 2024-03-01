# Mechanical oscillator in heat bath

using QuantumOptics, Plots, LaTeXStrings

# Parameters
ω_mech = 2π*1.1
Q = 1e7
γ0 = ω_mech / Q

# Constants
ħ = 1.05457182e-34 * 1e9^2 / 1e6
kB = 1.380649e-23 * 1e9^2 / 1e6^2
T = 1e-3
n = Int(round(1/(exp(ħ*ω_mech/(kB*T))-1)))

# Hilbert space
b_mech = FockBasis(50)

# Operators
b = destroy(b_mech)
bt = create(b_mech)
bt2 = dagger(b)

# Mechanical oscillator Hamiltonian
H_mech = ω_mech*bt*b


J = [b,bt]

rates = [γ0*(n+1),γ0*n]

# Initial conditions
Ψ0 = fockstate(b_mech,0)
tspan = collect(0:100:1e7)
dt = 0.1
@time tout, ρt = timeevolution.master(tspan,Ψ0, H_mech, J; dt = dt, rates = rates);

plot(tspan,real(expect(bt*b,ρt)),
    title = "Mechanical oscillator",
    xlabel = L"t [\mu s]",
    ylabel = L"\langle n_m \rangle",
    legend = false,
    lc = 2,
    dpi = 1400
)





