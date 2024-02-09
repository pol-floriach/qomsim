using DifferentialEquations, StaticArrays, Plots, Statistics, LaTeXStrings

# Energy of mechanical oscillator
function energy(sol)
    Q = real( @view sol[1,1:100:end])
    P = real(@view sol[2,1:100:end])
    return 1/2 * (P.^2 + Q.^2) 
end

# Harmonic oscillaotr (membrane)
function mean_posmom(u,p,t)
    @views Ω, Γ = p
    @views Q, P = u

    dQ = Ω*P
    dP = -Ω*Q - Γ*P
    SA[dQ, dP]
end

# ODE for light inside cavity only
function mean_light(α,p,t)
    @views κby2, iΔ, sqrtκinαin = p
    -(κby2 +iΔ)*α + sqrtκinαin
end

# Light + mechanical oscillator
function mean_posmomlight(u,p,t)
    @views Ω, Γ, κby2, iΔ, sqrt2g₀, sqrtκinαin = p
    @views Q, P, α = u

    dQ = Ω*P
    dP = -Ω*Q - Γ*P - sqrt2g₀*abs2(α)
    dα = -(κby2 + iΔ + im*sqrt2g₀*Q)*α + sqrtκinαin
    SA[dQ, dP, dα]
end

# Returns the mean α 
function alpha2simu(tspan, κ, κin, αin)
    Δvec = (-5e8:1e7:5e8) * im
    αss = zeros(ComplexF64,size(Δvec))

    for j in eachindex(Δvec)
        p = (κ/2, Δvec[j], sqrt(κin)*αin)
        prob = ODEProblem(mean_light!, 0.0im, tspan, p)
        sol_light = solve(prob,Tsit5(), saveat = 1e-8, dt = 1e-7, dtmin = 1e-10, dtmax = 1e-8, maxiters = tspan[2]/1e-10)
        αss[j] = mean(sol_light.u[Int(round(end/2)):end])
    end
    return αss
end


# Constants
begin
    Ω = 2π * 1.1e6
    Q_factor = 1e8 
    Γ = Ω / Q_factor

    κ = 2*pi*20e6
    Δ = 0.0
    iΔ = im*Δ

    g₀ = 2*pi*100
    κin = κ/100
    αin = 1e7
end
param = (Ω, Γ, κ/2, im*Δ, sqrt(2)*g₀, sqrt(κin)*αin)
#p = (κ/2, iΔ, sqrt(κin)*αin) # if only for 

u₀ = @SVector [1.0,0.0, 0.0im];

tspan = (0,0.01)

prob = ODEProblem(mean_posmomlight, u₀, tspan, param)
@time sol = solve(prob, Tsit5(), saveat = 1e-8, dt = 1e-10, dtmin = 1e-12, dtmax = 5e-10, maxiters = tspan[2]/1e-10)

N = Int64(round(tspan[2]/2/1e-6))
N = 1
tvec = sol.t[N:N+1000]
posvec = real(sol[1,N:N+1000])
momvec = real(sol[2,N:N+1000])
αvec = sol[3,N:N+1000]

p1 = plot(tvec,posvec, legend = false, xlabel = "t", ylabel = "x")#, xlims = (0,10/Ω))#ylims = (-1,1), xlims = ())
E = energy(sol)
p2 = plot(sol.t[1:1000:end], E[1:10:end], xlabel = "t", ylabel = "Energy", legend = false)

plot(p1,p2,layout = (2,1))

p3 = plot(tvec, abs2.(αvec), legend = false, xlabel = "t")


# What is the steady-state value?
function mean_posmomlight_steadystate(u,p)
    @views Ω, Γ, κby2, iΔ, sqrt2g₀, sqrtκinαin = p
    @views Q, α = u

    dP = -Ω*Q - sqrt2g₀*abs2(α)
    dα = -(κby2 + iΔ + im*sqrt2g₀*Q)*α + sqrtκinαin
    SA[dP, dα]
end
u₀_newton = SA[0.0, 0.0im]
prob_ss = NonlinearProblem(mean_posmomlight_steadystate, u₀_newton, param)
steady_state = solve(prob_ss, NewtonRaphson())






# # SOLVE ONLY LIGHT ODE
# prob = ODEProblem(mean_light, 0.0im, tspan, p)
# @time sol_light = solve(prob,Tsit5(), saveat = 1e-8, dt = 1e-7, dtmin = 1e-10, dtmax = 1e-8, maxiters = tspan[2]/1e-10)

# # SOLVING LIGHT ODE IN A RANGE OF Δ AND COMPARE mean(α) WITH STEADY-STATE VALUE
# α_simu = alpha2simu(tspan,κ, κin, αin)
# myrange = -5e8:1e4:5e8
# vals = @. αin * sqrt(κin) / (κ/2 + im * myrange)
# plot(myrange, abs2.(vals), label = "Steady-state", xlabel = L"\Delta", ylabel = L"|\alpha|^2", size = (650,400), dpi = 700)
# plot!(-5e8:1e7:5e8, abs2.(α_simu), lt = :scatter, markersize = 2, markerstrokewidth = 0, label = "Mean")
# savefig("N_vs_delta")


# SOLVE MEMBRANE EQUATIONS (damped harmonic)
# prob = ODEProblem(mean_posmom, u₀, tspan, p)
# @time sol = solve(prob, Tsit5(), saveat = 1e-6, dt = 1e-7, dtmin = 1e-10, dtmax = 1e-8, maxiters = tspan[2]/1e-10)
# Plot of a part of the time series (in half)
# N = Int64(round(tspan[2]/2/1e-6))
# tvec = sol.t[N:N+100]
# solvec = sol[1,N:N+100]
# p1 = plot(tvec,solvec, legend = false, xlabel = "t", ylabel = "x");
# # Energy calculation
# E = energy(sol)
# p2 = plot(sol.t[1:1000:end], E2[1:10:end], xlabel = "t", ylabel = "Energy", legend = false);
# plot(p1,p2, layout = (2,1))



