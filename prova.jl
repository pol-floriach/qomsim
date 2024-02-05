using DifferentialEquations, Plots, StaticArrays

function classical_membrane(u, p, t)
    Γ₀, Ω₀2, noiseconst = p
    @views x,z = u
    # x = u[1]
    # z = u[2]

    dx = z
    d2x = -Γ₀*z - Ω₀2*x
    SA[dx,d2x]
end

function classical_membrane!(du,u,p,t)
    @views Γ₀, Ω₀2, noiseconst = p
    @views x, z = u

    du[1] = z
    du[2] -Γ₀*z - Ω₀2*x
end

#Constants and initial conditions
begin
    Ω₀ = 2π*40e3
    const Q = 3e7;
    Γ₀ = Ω₀0 / Q

    const kB = 1.380649e-23;
    const m = 12e-9;
    const T = 5e-3;

    # Γ₀ = Γ₀0 / Ω₀0^2
    # Ω₀ = 1


    p = (Γ₀, Ω₀^2, sqrt(2*kB*T*m*Γ₀)/m);

    u₀ = @SVector rand(2)
    u₀ = u₀ / 100
    u₀2 = rand(2)
    # tspan = (0,1/Γ₀)
    tspan = (0,50)
end

# Without noise
prob = ODEProblem(classical_membrane, u₀, tspan, p);
@time sol = solve(prob, AutoVern9(Rodas5()), saveat = 1, dt = 1e-3, dtmin = 1e-7, dtmax = 1e-1,maxiters = tspan[2]/1e-15);

# Is it a stiff problem?
prob2 = ODEProblem(classical_membrane!,u₀2, tspan, p)
@time sol2 = solve(prob2, FBDF(), saveat = 1, dt = 1e-2, dtmin = 1e-5)


plot(sol.t, sol[1,:],
    xlabel = "t (s)",
    ylabel = "x",
    legend = false,
    title = "Classical SHO time series")



diff= @. max(abs(sol[1,:]- sol2[1,:]))

