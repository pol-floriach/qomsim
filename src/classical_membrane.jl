using DifferentialEquations, Plots, StaticArrays, LinearAlgebra, ProgressLogging, Statistics

# function classical_membrane(u, p, t)
#     @views @inbounds Γ₀, Ω₀2 = p[1:2]
#     @views x,v = u
#     dx = v
#     d2x = - Ω₀2*x  -Γ₀*v
#     SA[dx,d2x]
# end

function classical_membrane(du, u, p, t)
    @views @inbounds Γ₀, Ω₀2 = p[1:2]
    @views x,v = u

    du[1] = v
    du[2] = - Ω₀2*x  -Γ₀*v
    nothing
end


# function nondimensionalized(u,p,t)
#     @views γ = p[1]
#     @views x,v = u

#     dx = v
#     d2x = -γ*v - x

#     SA[dx,d2x]
# end

# function noise_fun(u,p,t)
#     return [0.0;@views p[3]]
# end

function noise_fun(du,u,p,t)
    du[1] = 0.0
    du[2] = p[3]
end


meanx2(T) = 2*ħ*Ω₀/(kB*T)*sqrt(ħ/(2*m*Ω₀))


#Constants and initial conditions
begin
    Ω₀ = 2π * 1.1e6
    Q = 1e8 
    Γ₀ = Ω₀ / Q

    const kB = 1.380649e-23;
    const m = 12e-9;
    T = 273#5e-3;
    const ħ = 1.05457182e-34

    p = (Γ₀, Ω₀^2, sqrt(2*kB*T*m*Γ₀/xzp)/m);
    #p = (Γ₀, Ω₀^2,0.0);
    p_nondim = (1/Q,1/sqrt(Ω₀));

    u₀ =  [1,0];
    u₀s = SA[u₀...];
   
    u₀_nondim = SA[0, m*sqrt(Ω₀)/sqrt(2*kB*T*m*Γ₀)]
    #u₀_nondim = SA[m*Ω₀*sqrt(Ω₀)/sqrt(2*kB*T*m*Γ₀),0]
end



tspan = (0,1e-3)
# Without noise (and using StaticArrays)
prob = ODEProblem(classical_membrane, u₀, tspan, p)
@time sol = solve(prob, Tsit5(), saveat = 1e-8, dt = 1e-9, dtmin = 1e-10, dtmax = 1e-8, maxiters = tspan[2]/1e-12, progress = true);

# Now with NOISE term (Fo * ξ(t)~Normal(0,1))
probnoise = SDEProblem(classical_membrane, noise_fun, u₀, tspan,p)
@time sol_stoch = solve(probnoise,EulerHeun(), dt = 1e-10, dtmax = 1e-8, saveat = 5e-9, reltol = 1e-8, abstol = 1e-8,progress = true, maxiters = tspan[2]/1e-12)


N = size(sol,2)
N2 = size(sol_stoch,2)
# myrange = Int(round(N/2))-100000:Int(round(N/2))
myrange =  1:100000
# myrange2 = 1:100000
myrange = N-1000:N
myrange2 = N2-2000:N2
plot(sol.t[myrange], sol[1,myrange], xlabel = "t (s)", ylabel = "x", label = "No Noise", title = "Classical SHO time series", ls = :dash)
plot!(sol_stoch.t[myrange2],sol_stoch[1,myrange2], label = "Noise")

#plot(sol.t[1:1000:end],abs.(sol[1,1:1000:end] .- sol_stoch[1,1:10000:end]))

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
prob_nondim_stoch = SDEProblem(nondimensionalized, noise, u₀_nondim, τspan, p_nondim)
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

