using DifferentialEquations, Plots, StaticArrays

function classical_membrane(u, p, t)
    Γ₀, Ω₀2, noiseconst = p
    x = u[1]
    z = u[2]

    dx = z
    d2x = -Γ₀*z - Ω₀2 *x
    SA[dx,d2x]
end


function classical_feedback(du,u,p,t)
    Γ₀, Ω₀2, noiseconst = p
    x = u[1]
    z = u[2]
    y = u[3]
    v = u[4]

    du[1] = z
    du[2] = -Γ₀*du[1] - Ω₀2*x - g*Γ₀*v
    du[3] = v
    nothing
end


begin
    const Ω₀ = 2π*40e3
    const Q = 3e7;
    const Γ₀ = Ω₀ / Q

    const kB = 1.380649e-23;
    const m = 12e-9;
    const T = 5e-3;

    p = (Γ₀, Ω₀^2, sqrt(2*kB*T*m*Γ₀)/m);

    u₀ = @SVector rand(2)
    u₀ = u₀ / 100
    # tspan = (0,1/Γ₀)
    tspan = (0,100)
end

# Without noise
prob = ODEProblem(classical_membrane, u₀, tspan, p);

@time sol = solve(prob, Vern9(), dt = 1e-10, dtmax = 1e-3, dtmin = 1e-16, maxiters = tspan[2]/1e-15, saveat =  1e-1)
# @time sol = solve(prob, Tsit5(), dt = 1e-6, dtmax = 1e-3, dtmin = 1e-15, maxiters = tspan[2]/1e-15, saveat =  1e-2);

plot(sol.t, sol[1,:])


# With noise
function noise_fun(u,p,t)
    return [0;p[3]]
end;
prob2 = SDEProblem(classical_membrane, noise_fun, u₀, tspan,p);
sol_stoch = solve(prob2, dt =1e-10, dtmax = 1e-3, dtmin = 1e-16, maxiters = tspan[2]/1e-15, saveat = 1e-2, progress = true)
plot(sol_stoch.t,sol_stoch[1,:]);


# Stability?
using LinearAlgebra
J = [0 1; -Ω₀^2 -Γ₀]
eigvals(J)


# Comparison with theoretical
using FFTW, LinearAlgebra
dft_sol = abs.(fft(sol[1,:]))

Ω = collect(1:length(dft_sol))
plot(Ω, dft_sol)


const xzp = 4e-15
const nth = kB*T/(m*Γ₀)
const ħ = 1.054571817e-34
g = 1.4e5

χg = @. 1 / ((1+g) + 2im*(Ω - Ω₀)/Γ₀)
Sxxzp = 4 * xzp^2 / Γ₀
#Sxximp = 2* ħ^2 * Q / (kB*T*m*Ω₀)
Sxximp = 1e-17

nimp = Sxximp / (2*Sxxzp)

Sxx = @. 2 * Sxxzp * (nth + g^2*nimp) * abs(χg)^2

plot(Sxx)