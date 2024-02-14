# Script for running semiclassical equations (2.28) in a range of Δ
using DifferentialEquations, StaticArrays, Plots, Statistics, ProgressLogging

# Light + mechanical oscillator
function mean_posmomlight(u,p,t)
    @views Ω, Γ, κby2, iΔ, sqrt2g₀, sqrtκinαin = p
    @views Q, P, α = u

    dQ = Ω*P
    dP = -Ω*Q - Γ*P - sqrt2g₀*abs2(α)
    dα = -(κby2 + iΔ + im*sqrt2g₀*Q)*α + sqrtκinαin
    SA[dQ, dP, dα]
end

function α2forΔrange(u₀, tspan,Δt; Δvec = (-3e8:1e7:3e8))
    # Constants
    begin
        Ω = 2π * 1.1e6
        Q_factor = 1e8 
        Γ = Ω / Q_factor
    
        κ = 2*pi*20e6
    
        g₀ = 2*pi*100
        κin = κ/100
        αin = 1e7
    end;

    iΔvec = im*Δvec
    αmat = zeros(ComplexF64, size(Δvec,1),Int(tspan[2]/Δt)+1)
    
    for j in eachindex(Δvec)
        param = (Ω, Γ, κ/2, im*Δvec[j], sqrt(2)*g₀, sqrt(κin)*αin)
        prob = ODEProblem(mean_posmomlight, u₀, tspan, param)
        sol = solve(prob, Tsit5(), saveat = Δt, dt = 1e-10, dtmin = 1e-11, dtmax = 1e-8, maxiters = tspan[2]/1e-10)
        println(j)
        αmat[j,:] = sol[1,:]
        # println(size(sol[3,:]))
    end


    return αmat
end


u₀ = @SVector [0.0,0.0, 0.0im];
tspan = (0,1e-5)
Δt = 1e-9
Δvec = (-2e8:5e6:2e8)
myαmat = α2forΔrange(u₀, tspan, Δt;  Δvec = Δvec);
timevec = collect(tspan[1]:Δt:tspan[2])

myrange = 1:100
heatmap(timevec[myrange], 
        Δvec,
        (abs2.(myαmat[:,myrange])),
        xlabel = "t",
        ylabel = "Δ",
        colormap = :viridis,
        )

idxs = collect(1:40)
idxs = [1,10,41,70,80]
myrange = 1:3000

pp = plot();
for i in idxs
    plot!(timevec[myrange], real.(myαmat[i,myrange]), label = "Δ = $(Δvec[i])", legend = :topright) #lt = :scatter, markersize = 2, markerstrokewidth = 0)#(19488.5, 19491))
end
pp