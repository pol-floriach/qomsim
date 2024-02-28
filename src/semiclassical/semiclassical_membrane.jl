using DifferentialEquations, StaticArrays, Plots, Statistics, ProgressLogging, Peaks

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
    dP = -Ω*Q -Γ*P
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


# Returns mean(α) for different Δ values 
# function alpha2simu(tspan, κ, κin, αin)
#     Δvec = (-5e8:1e7:5e8) * im
#     αss = zeros(ComplexF64,size(Δvec))

#     for j in eachindex(Δvec)
#         p = (κ/2, Δvec[j], sqrt(κin)*αin)
#         prob = ODEProblem(mean_light!, 0.0im, tspan, p)
#         sol_light = solve(prob,Tsit5(), saveat = 1e-8, dt = 1e-7, dtmin = 1e-10, dtmax = 1e-8, maxiters = tspan[2]/1e-10)
#         αss[j] = mean(sol_light.u[Int(round(end/2)):end])
#     end
#     return αss
# end


# Constants
begin
    Ω = 2π * 1.1e6
    Q_factor = 1e8 
    Γ = Ω / Q_factor

    κ = 2*pi*20e6
    Δ = 0.5e8
    iΔ = im*Δ

    g₀ = 2*pi*100
    κin = κ/100
    αin = 1e7
end
p = (Ω,Γ)
param = (Ω, Γ, κ/2, im*Δ, sqrt(2)*g₀, sqrt(κin)*αin)
#p = (κ/2, iΔ, sqrt(κin)*αin) # if only for light

u₀ = @SVector [1.0,0.0, 0.0im];
# u₀ = @SVector [1.0, 0.0]

tspan = (0.0,1e-5)

prob = ODEProblem(mean_posmomlight, u₀, tspan, param)
@time sol = solve(prob, Tsit5(), saveat = 1e-9, dt = 1e-10, dtmin = 1e-11, dtmax = 1e-8, maxiters = tspan[2]/1e-12, progress = true)

pksmaxq, valsmaxq = findmaxima(real.(sol[1,400:end-1000]));
pksmaxa, valsmaxa = findmaxima(real(sol[3,400:end-1000]));

plot_q = plot(sol.t, real(sol[1,400:end-1000]), lt = :scatter, markersize = 1, markerstrokewidth = 0)
plot!(sol.t[pksmaxq],valsmaxq, lt = :scatter)

plot_a = plot(sol.t, real.(sol[3,400:end-1000]), lt = :scatter, markersize = 1, markerstrokewidth = 0)#(19488.5, 19491))
plot!(sol.t[pksmaxa],valsmaxa, lt = :scatter)

plot(plot_q, plot_a, layout = (2,1))

freq_q = 0.0
for i in 2:size(pksmax,1)
freq_q += pksmaxq[i] - pksmaxq[i-1]
end
freq_q = 1/(freq_q / (size(pksmaxq,1)-1) * 1e-9)

freq_a = 0.0
for i in 2:size(pksmax,1)
freq_a += pksmaxa[i] - pksmaxa[i-1]
end
freq_a = 1/(freq_a / (size(pksmaxa,1)-1) * 1e-9)



N = Int64(round(tspan[2]/1e-7))
# myrange = N-100:N
myrange = 1:100
αrange = 1:10
tvec = sol.t[myrange]
posvec = real(sol[1,myrange])
momvec = real(sol[2,myrange])
αvec = sol[3,αrange]

p1 = plot(tvec,posvec, legend = false, xlabel = "t", ylabel = "x", xlims = (0,2π*5/Ω))
p2 = plot(tvec,momvec, legend = false, xlabel = "t", ylabel = "p", xlims = (0,2π*5/Ω))

p3 = scatter(sol.t[αrange], abs2.(αvec), legend = false, xlabel = "t", ylabel = "|α|²")

E2 = energy(sol)
p4 = plot(sol.t[1:10000:end], E2[1:100:end], xlabel = "t", ylabel = "Energy", label = "E", xlims = (0,0.5))

plot!(sol.t[1:1000:end], abs2.(sol[1,1:1000:end])/2, label = "Q²")
plot!(sol.t[1:1000:end], abs2.(sol[2,1:1000:end])/2, label = "P²")



# Maybe find maxima?
pksmax, valsmax = findmaxima(real(sol[1,1:10:end]));
pksmin, valsmin = findminima(real(sol[1,1:10:end]));
p5 = plot(sol.t[pksmax[1:1000:end]*10], valsmax[1:1000:end], lw = 2, xlabel = "t", ylabel = "x")
plot!(sol.t[pksmin[1:1000:end]*10], valsmin[1:1000:end], fillrange = valsmax[1:1000:end], fillalpha = 0.35, c = 3, legend = false, lw = 2)
plot!(sol.t[pksmin[1:1000:end]*10], zeros(size(pksmin[1:1000:end])), c = :black, lw = 1.5)


pksmaxp, valsmaxp = findmaxima(real(sol[2,1:10:end]));
pksminp, valsminp = findminima(real(sol[2,1:10:end]));
p6 = plot(sol.t[pksmaxp[1:1000:end]*10], valsmaxp[1:1000:end], lw = 2, xlabel = "t", ylabel = "p")
plot!(sol.t[pksminp[1:1000:end]*10], valsminp[1:1000:end], fillrange = valsmaxp[1:1000:end], fillalpha = 0.35, c = 2, legend = false, lw = 2)
plot!(sol.t[pksminp[1:1000:end]*10], zeros(size(pksminp[1:1000:end])), c = :black, lw = 1.5)



plot(p1,p4, p2, p5, p3, p6, layout = (3,2), plot_title = "Δ = $(Δ)", size = (900,600), dpi = 800)

plot(p4,p5, p6, layout = (3,1), plot_title = "Δ = $(Δ)", size = (700,550), dpi = 800)


codi
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



