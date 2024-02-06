using DifferentialEquations, Plots, StaticArrays, LinearAlgebra

function classical_membrane(u, p, t)
    @views Γ₀, Ω₀2 = p[1:2]
    @views x,v = u
    dx = v
    d2x = -Γ₀*v - Ω₀2*x
    SA[dx,d2x]
end

function nondimensionalized(u,p,t)
    @views γ = p[1]
    @views x,v = u

    dx = v
    d2x = -γ*v - x

    SA[dx,d2x]
end

# function classical_membrane!(du,u,p,t)
#     @views Γ₀, Ω₀2 = p[1:2]
#     @views x, v = u

#     du[1] = v
#     du[2] = -Γ₀*v - Ω₀2*x
# end

function noise_fun(u,p,t)
    return [0; @views p[2]]
end

function noise_nondim(u,p,t)
    return[0;1]
end

#Constants and initial conditions
begin
    Ω₀ = 2π;
    Q = 1e3;
    Γ₀ = Ω₀ / Q
    

    const kB = 1.380649e-23;
    const m = 12e-9;
    const T = 5e-3;
   

    p = (Γ₀, Ω₀^2, sqrt(2*kB*T*m*Γ₀)/m);
    p_nondim = (1/Q,sqrt(2*kB*T*m*Γ₀)/m);

    u₀ =  [1,0];
    u₀s = SA[u₀...];
   
    u₀_nondim = SA[0, m*Ω₀/sqrt(2*kB*T*m*Γ₀)]

    tspan = (0,1/Γ₀);
    τspan = (0,Ω₀/Γ₀);
end


# Without noise (and using StaticArrays)
prob = ODEProblem(classical_membrane, u₀s, tspan, p);
@time sol = solve(prob, AutoTsit5(Rodas5()), saveat = 0.01, dt = 1e-3, dtmin = 1e-7, dtmax = 1e-2);


# Now with NOISE term (Fo * ξ(t)~Normal(0,1))
probnoise = SDEProblem(classical_membrane, noise_fun, u₀s, tspan,p);
@time sol_stoch = solve(probnoise, dt =1e-5, dtmax = 1e-2, dtmin = 1e-7, saveat = 0.01);#, maxiters = tspan[2]/1e-15, progress = true)


# Analytical solution
dt = 0.01;
begin
    N = 10000
    Γ = Γ₀
    Ω = sqrt(Ω₀^2 - Γ₀^2)
    t = dt*collect(0:N-1)
    xtheory = @. exp(-0.5*Γ*t) * sin(Ω*t)

    λ₊ = -Γ/2 + im*Ω
    λ₋ = conj(λ₊)

    # THE MAGICAL CORRECTIONS
    Γprime = -2*real(exp(λ₊*dt)-1)/dt 
    Ωprime = imag(exp(λ₊*dt)-1)/dt
    Φ = [1-0.5*Γprime*dt Ωprime*dt ; -Ωprime*dt 1-0.5*Γprime*dt]


    xs = zeros(N,2)
    xs[1,:] = [1;0]
    for i in 2:N
        xs[i,:] = Φ*xs[i-1,:]
    end
end

 
# Non-dimensionalized version 
prob_nondim = ODEProblem(nondimensionalized, u₀_nondim, τspan, p_nondim)
@time sol_nondim = solve(prob_nondim, saveat = 0.01*Ω₀, dt = 1e-3, dtmin = 1e-7, dtmax = 1e-2)
# Non-dimensionalized with noise
prob_nondim_stoch = SDEProblem(nondimensionalized, noise_nondim, u₀_nondim, τspan, p_nondim)
@time sol_nondim_stoch = solve(prob_nondim_stoch, saveat = 0.01*Ω₀ , dt = 1e-3)#, dtmin = 1e-7, dtmax = 1e-2)



# TIME-SERIES PLOT
plot(sol.t, sol[1,:], xlabel = "t (s)", ylabel = "x", label = "No Noise", title = "Classical SHO time series", xlims = (0,5), ylims = (-1.25,1.25))

plot!(sol_nondim.t/Ω₀,sol_nondim[1,:]*sqrt(2*kB*T*m*Γ₀)/(m*Ω₀^2), xlims = (0,5), label ="Non-dimensonalized (no noise)", ls = :dash)

plot!(t,xs[:,1], label = "Analytical", xlims = (0,5))

plot!(sol_stoch.t,sol_stoch[1,:], label = "w/ Noise")

plot!(sol_nondim_stoch.t/Ω₀,sol_nondim_stoch[1,:]*sqrt(2*kB*T*m*Γ₀)/(m*Ω₀^2), xlims = (0,5), label ="Non-dimensonalized (w/ noise)", ls = :dash)



# Error from dimensional and non-dimensional simulation
max_abs_err = max(@. abs(sol[1,:] - sol_nondim[1,:]*sqrt(2*kB*T*m*Γ₀)/(m*Ω₀^2))...)

# Energy conservation? 

function energy(sol)
    @. E = 1/2*m*sol.u[1,:]^2 
end