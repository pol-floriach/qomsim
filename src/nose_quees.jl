using DifferentialEquations, Plots, StaticArrays, LinearAlgebra, ProgressLogging, Statistics

function classical_membrane(u, p, t)
    @views Γ₀, Ω₀2 = p[1:2]
    @views x,v = u
    dx = v
    d2x = - Ω₀2*x  -Γ₀*v
    SA[dx,d2x]
end


function thermal_noise(u,p,t)
    SA[0.0;@views p[3]]
end

function back_action(u,p,t)
    SA[0.0; @views p[4]]
end


#Constants (note that we have chosen [x] = nm and [t] = μs!)
begin
    Ω₀ = 2π*1.1                         # [Ω₀] = s^(-1)
    Q = 1e7
    Γ₀ = Ω₀ / Q                         # [Γ₀] = s^(-1)

    kB = 1.380649e-23  * 1e9^2/1e6^2    # [kB] = nm^2 kg  μs^(-2) K^(-1)
    m = 1e-12;                          # [m]  = kg
    T = 0;                            # [T]  = K
    ħ = 1.05457182e-34  * 1e9^2 / 1e6   # [ħ]  = m^2 kg / s

    xzpf = sqrt(ħ/(2*m*Ω₀))           
    nth = 1/(exp(ħ*Ω₀/(kB*T)) - 1)
end




# Parameters and initial conditions
begin
    p = (Γ₀, Ω₀^2, sqrt(2*kB*T*m*Γ₀)/m, 2*ħ*2π/λ*sqrt(Nflux));
    u₀ =  [xzpf * sqrt(2*nth),0];
    u₀s = SA[u₀...];
end



tspan = (0,1e3)
# Without noise (and using StaticArrays)
prob = ODEProblem(classical_membrane, u₀s, tspan, p)
@time sol = solve(prob, Tsit5(), saveat = 1e-1, dt = 1e-3, dtmin = 1e-7, dtmax = 1e-2, maxiters = tspan[2]/1e-12, progress = true, reltol = 1e-10, abstol = 1e-10);

# Now with NOISE term (Fo * ξ(t)~Normal(0,1)) 
probnoise = SDEProblem(classical_membrane, back_action, u₀s, tspan,p)
@time sol_stoch = solve(probnoise,dt = 1e-4, dtmax = 1e-4, saveat = 1e-1, progress = true, maxiters = tspan[2]/1e-8)

plot(sol.t, sol[2,:], xlabel = "t (μs)", ylabel = "x (nm)", label = "w/o noise", ls = :dash);
plot!(sol_stoch.t, sol_stoch[2,:], label = "w/ Noise", ls = :dashdot)


println("Deterministic σₓ = $(std(sol[1,Int(round(end/2)):end])) nm")
println("Stochastic σₓ = $(std(sol_stoch[1,Int(round(end/2)):end])) nm")
println("Theoretical σₓ = $(sqrt(2*nth)*xzpf) nm")

N = size(sol,2)
N2 = size(sol_stoch,2)
# myrange = Int(round(N/2))-100000:Int(round(N/2))

# myrange = N-100000:N
# myrange2 = N2-100000:N2
myrange = 1:100:N
myrange2 = 1:100:N2
plot(sol.t[myrange], sol[2,myrange], xlabel = "t (μs)", ylabel = "x(nm)", label = "No Noise", title = "Classical SHO time series", ls = :dash)
plot!(sol_stoch.t[myrange2],sol_stoch[2,myrange2], label = "Noise")

p2 = plot(sol_stoch.W[2,1:100:end]*p[3])

plot(p1,p2, layout = (2,1) )




# NOT USED AS OF RIGHT NOW

# # Analytical solution
# dt = 1e-4;
# begin
#     N = 1000000
#     Γ = Γ₀
#     Ω = sqrt(Ω₀^2 - Γ₀^2)
#     t = dt*collect(0:N-1)
#     xtheory = @. exp(-0.5*Γ*t) * sin(Ω*t)

#     λ₊ = -Γ/2 + im*Ω
#     λ₋ = conj(λ₊)

#     # THE MAGICAL CORRECTIONS
#     Γprime = -2*real(exp(λ₊*dt)-1)/dt 
#     Ωprime = imag(exp(λ₊*dt)-1)/dt
#     Φ = [1-0.5*Γprime*dt Ωprime*dt ; -Ωprime*dt 1-0.5*Γprime*dt]
#     Φ = [1 dt; -Ωprime^2*dt 1-Γprime*dt]

#     xs = zeros(N,2)
#     xs[1,:] = [0;1]
#     for i in 2:N
#         xs[i,:] = Φ*xs[i-1,:]
#     end
# end

# plot(t,xs[:,1], label = "Analytical")
 

# # Non-dimensionalized version 
# prob_nondim = ODEProblem(nondimensionalized, u₀_nondim, τspan, p_nondim)
# @time sol_nondim = solve(prob_nondim,Vern9(), saveat = 0.01*Ω₀, dt = 1e-3, dtmin = 1e-5, dtmax = 1e-1)
# # Non-dimensionalized with noise
# prob_nondim_stoch = SDEProblem(nondimensionalized, noise, u₀_nondim, τspan, p_nondim)
# @time sol_nondim_stoch = solve(prob_nondim_stoch, EM(), saveat = 0.01*Ω₀, dt = 1e-3)#, dtmin = 1e-7, dtmax = 1e-2)



# # TIME-SERIES PLOT
# plot(sol.t, sol[1,:], xlabel = "t (s)", ylabel = "x", label = "No Noise", title = "Classical SHO time series")#, ylims = (-1,1))

# plot!(sol_nondim.t/Ω₀,sol_nondim[1,:]*sqrt(2*kB*T*m*Γ₀)/(m*Ω₀^(3/2)), label ="Non-dimensonal (no noise)", ls = :dash)

# plot!(t,xs[:,1], label = "Analytical")

# plot!(sol_stoch.t,sol_stoch[1,:], label = "w/ Noise")

# plot!(sol_nondim_stoch.t/Ω₀,sol_nondim_stoch[1,:]*sqrt(2*kB*T*m*Γ₀)/(m*Ω₀*sqrt(Ω₀)), label ="Non-dimensonal (w/ noise)", ls = :dash)



# # Error from dimensional and non-dimensional simulation
# max_abs_err = maximum(@. abs(sol[1,:] - sol_nondim[1,:]*sqrt(2*kB*T*m*Γ₀)/(m*Ω₀^2)))
# max_abs_err = maximum(@. abs(sol[1,1:N] - xs[:,1]))

# # Energy conservation? 

# function energy(sol)
#      E = @. 1/2*m*sol[1,:]^2 
# end

# plot(sol.t,energy(sol_nondim))

