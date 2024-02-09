using DifferentialEquations, Plots, StaticArrays, LinearAlgebra

function classical_membrane(u, p, t)
    @views Γ₀, Ω₀2 = p[1:2]
    @views x,v = u
    dx = v
    d2x = - Ω₀2*x  -Γ₀*v
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
    # return [0; @views p[2]]
    return [0.0;p[2]]
end

function noise_nondim(u,p,t)
    # return[0; @views p[2]]
    return [0.0; p[2]]
end

#Constants and initial conditions
begin
    Ω₀ = 2π*40e3;
    Q = 1e7;
    Γ₀ = Ω₀ / Q
    

    const kB = 1.380649e-23;
    const m = 12e-9;
    const T = 5e-3;
   

    p = (Γ₀, Ω₀^2, sqrt(2*kB*T*m*Γ₀)/m);
    p_nondim = (1/Q,1/sqrt(Ω₀));

    u₀ =  [1,0];
    u₀s = SA[u₀...];
   
    u₀_nondim = SA[0, m*sqrt(Ω₀)/sqrt(2*kB*T*m*Γ₀)]
    #u₀_nondim = SA[m*Ω₀*sqrt(Ω₀)/sqrt(2*kB*T*m*Γ₀),0]
    tspan = (0,1/Γ₀);
    τspan = (0,Ω₀/Γ₀);
end

tspan = (0,10)
# Without noise (and using StaticArrays)
prob = ODEProblem(classical_membrane, u₀s, tspan, p);
@time sol = solve(prob, Tsit5(), saveat = 1e-6, dt = 1e-3, dtmin = 1e-10, dtmax = 1e-6, maxiters = tspan[2]/1e-15);
#@time sol = solve(prob, VerletLeapfrog(), saveat = 0.01, dt = 1e-3, dtmin = 1e-7, dtmax = 1e-2, maxiters = tspan[2]/1e-15);
plot(sol.t, sol[1,:], xlabel = "t (s)", ylabel = "x", label = "No Noise", title = "Classical SHO time series", xlims = (0,100e-5), ylims = (-1,1))#, ylims = (-1,1))



# Now with NOISE term (Fo * ξ(t)~Normal(0,1))
probnoise = SDEProblem(classical_membrane, noise_fun, u₀s, tspan,p);
@time sol_stoch = solve(probnoise, EM(), dt =1e-5, dtmax = 1e-2, dtmin = 1e-7, saveat = 0.01)
#,


# Analytical solution
dt = 1e-4;
begin
    N = 1000000
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
    Φ = [1 dt; -Ωprime^2*dt 1-Γprime*dt]

    xs = zeros(N,2)
    xs[1,:] = [0;1]
    for i in 2:N
        xs[i,:] = Φ*xs[i-1,:]
    end
end

plot(t,xs[:,1], label = "Analytical")
 

# Non-dimensionalized version 
prob_nondim = ODEProblem(nondimensionalized, u₀_nondim, τspan, p_nondim)
@time sol_nondim = solve(prob_nondim,Vern9(), saveat = 0.01*Ω₀, dt = 1e-3, dtmin = 1e-5, dtmax = 1e-1)
# Non-dimensionalized with noise
prob_nondim_stoch = SDEProblem(nondimensionalized, noise_nondim, u₀_nondim, τspan, p_nondim)
@time sol_nondim_stoch = solve(prob_nondim_stoch, EM(), saveat = 0.01*Ω₀, dt = 1e-3)#, dtmin = 1e-7, dtmax = 1e-2)



# TIME-SERIES PLOT
plot(sol.t, sol[1,:], xlabel = "t (s)", ylabel = "x", label = "No Noise", title = "Classical SHO time series")#, ylims = (-1,1))

plot!(sol_nondim.t/Ω₀,sol_nondim[1,:]*sqrt(2*kB*T*m*Γ₀)/(m*Ω₀^(3/2)), label ="Non-dimensonal (no noise)", ls = :dash)

plot!(t,xs[:,1], label = "Analytical")

plot!(sol_stoch.t,sol_stoch[1,:], label = "w/ Noise")

plot!(sol_nondim_stoch.t/Ω₀,sol_nondim_stoch[1,:]*sqrt(2*kB*T*m*Γ₀)/(m*Ω₀*sqrt(Ω₀)), label ="Non-dimensonal (w/ noise)", ls = :dash)



# Error from dimensional and non-dimensional simulation
max_abs_err = maximum(@. abs(sol[1,:] - sol_nondim[1,:]*sqrt(2*kB*T*m*Γ₀)/(m*Ω₀^2)))
max_abs_err = maximum(@. abs(sol[1,1:N] - xs[:,1]))

# Energy conservation? 

function energy(sol)
     E = @. 1/2*m*sol[1,:]^2 
end

plot(sol.t,energy(sol_nondim))

