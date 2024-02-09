function mean_posmomlight_steadystate(u,p)
    @views Ω, Γ, κby2, iΔ, sqrt2g₀, sqrtκinαin = p
    @views Q, α = u

    dP = -Ω*Q - sqrt2g₀*abs2(α)
    dα = -(κby2 + iΔ + im*sqrt2g₀*Q)*α + sqrtκinαin
    SA[dP, dα]
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


u₀_newton = SA[0.0, 0.0im]
prob_ss = NonlinearProblem(mean_posmomlight_steadystate, u₀_newton, param)
steady_state = solve(prob_ss, NewtonRaphson(), reltol = 1e-10, abstol = 1e-10)

steady_state.stats

steady_state

u₀ = SA[real(steady_state[1]), 0.0, steady_state[2]]


# Check if solution makes sense:
qs_error =abs(-sqrt(2)*g₀*abs2(steady_state[2])/Ω -steady_state[1])
# It does, as substituting numerical |α_ss|^2 into analytical q_ss yields the same result as numerical