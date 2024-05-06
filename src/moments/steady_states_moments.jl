using DifferentialEquations, StaticArrays, Plots, LaTeXStrings

function covariances_ss(u,p)
    @views γ0, ω, k, η = p
    @views Vx, Vp, Cxp = u
    dVx   = -γ0*Vx + 2*ω*Cxp - 8*k*η*Vx^2  + 0.5*γ0   
    dVp   = -γ0*Vp - 2*ω*Cxp - 8*k*η*Cxp^2 + 0.5*γ0 + 2*k  
    dCxp  = ω*Vp - ω*Vx - γ0*Cxp - 8*k*η*Vx*Cxp
    SA[dVx, dVp, dCxp]
end
# Parameters
ω = 2π*1.1
Q = 1e7
γ0 = ω / Q
Vx_th = Vp_th = 0.5 # at T = 0K, n̄ = 0
η = 0.9


kvec = collect(0.0:1e-2:10)
steady_values = zeros(length(kvec), 3)
ii = 0
u0_newt = SA[Vx_th, Vp_th, 0.0]
for k in kvec
    ii+=1
    p = (γ0, ω, k, η)
    prob_ss = NonlinearProblem(covariances_ss, u0_newt, p)
    ss = solve(prob_ss, NewtonRaphson(), reltol = 1e-10, abstol = 1e-10)
    steady_values[ii,:] = ss.u
    u0_newt = ss.u
end

plot(kvec, steady_values[:,1], 
    label = L"V_x",
    xlabel = "k", 
    ylabel = "Covariance (unitless)",
    size = (800,400),
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
hline!([0.5], label = L"\frac{1}{2}", legend = :outertopright, ticks = :native)

# plot!(kvec, steady_values[:,3],
#     label = L"C_{xp}",
#     lt = :scatter,
#     ms = 1.5,
#     msw = 0,
#     ma = 0.5)
println("Valor màxim de V_x és: ", maximum(steady_values[:,1]), ", valor mínim: ", minimum(steady_values[:,1]))
println("Valor màxim de V_p és: ", maximum(steady_values[:,2]), ", valor mínim: ", minimum(steady_values[:,2]))

plot(ke