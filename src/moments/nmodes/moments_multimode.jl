# Script for time evolution of the moments and covariances (SDE system) including 3 modes.
# Note: X and Y are the slowly-varying operators of x and p
using DifferentialEquations, Plots, LaTeXStrings, ProgressLogging

# Function to obtain O_k index of matrix (in a row or column)
function right_index(O,k,target= 1)
    return O == target ?  2k-1 : 2k
end
# Function to get index of a matrix stored in an array in column-major order (for 6x6 matrix)
colmaj(row, col) = @. row + (col-1)*6

# Kronecker delta function for better legibility
δ(x,y) = ==(x,y)

# Time evolution w/o information gain. TODO If too many allocations: try with staticarray / add Mi, Nj etc as parameters and mutate
function moments_evolution_3modes(du, u, p, t)
    # Parameters
    @views γ0, ωm, ωs, Γba ,Γmeas = p
    # X, Y, C
    Xs = @view u[1:2:5]
    Ys = @view u[2:2:6]
    C = @view u[7:end]
    # dX, dY, dC
    dXs = @view du[1:2:5]
    dYs = @view du[2:2:6]
    dC = @view du[7:end]

    #First moments
    @. dXs = -γ0/2*Xs + (ωs - ωm)*Ys
    @. dYs = -γ0/2*Ys - (ωs - ωm)*Xs

    # Covariance matrix
    for M = 1:2, N = 1:2
        for i in 1:3, j in 1:3
            # Getting the right indices
            Mi = right_index(M,i)
            Nj = right_index(N,j)
            Mconji = right_index(M,i,2)
            Nconjj = right_index(N,j,2)

            dC[colmaj(Mi, Nj)] = - (γ0[i]+γ0[j])/2 * C[colmaj(Mi,Nj)] + δ(Mi,Nj)*γ0[i]*(nth+0.5) + δ(M,N)*sqrt(Γba[i]*Γba[j]) +
                     + (-1)^δ(M,2)*(ωs[i] - ωm)*C[colmaj(Mconji, Nj)] + (-1)^δ(N,2)*(ωs[j] - ωm)*C[colmaj(Mi, Nconjj)] +
                     -4*(sum([sqrt(Γmeas[k])*C[colmaj(Mi, 2k-1)] for k in 1:3]))*(sum([sqrt(Γmeas[l])*C[colmaj(Nj, 2l-1)] for l in 1:3])) +
                     -4*(sum([sqrt(Γmeas[k])*C[colmaj(Mi, 2k)] for k in 1:3]))*(sum([sqrt(Γmeas[l])*C[colmaj(Nj, 2l)] for l in 1:3]))
        end
    end
    nothing
end

function infogain_3modes(du,u,p,t)
    @views Γmeas = p[end]
    C = @view u[7:end]
    du[1:2:5,1] = 2*[sum([sqrt(Γmeas[j])*C[colmaj(Xi,2j-1)] for j in 1:3]) for Xi in 1:2:5] # dWx in d<Xi>
    du[1:2:5,2] = 2*[sum([sqrt(Γmeas[j])*C[colmaj(Xi,2j)]       for j in 1:3]) for Xi in 1:2:5] # dWy in d<Xi>
    du[2:2:6,1] = 2*[sum([sqrt(Γmeas[j])*C[colmaj(Yi,2j-1)] for j in 1:3]) for Yi in 2:2:6] # dWx in d<Yi>
    du[2:2:6,2] = 2*[sum([sqrt(Γmeas[j])*C[colmaj(Yi,2j)]       for j in 1:3]) for Yi in 2:2:6] # dWy in d<Yi>

    du[7:end,1] .= 0.0
    du[7:end,2] .= 0.0
end

# Parameters
m = 1e-12
ωm = 2π*1.1
ωs = [2π*1.1, 2π*2, 2π*0.2]
Qs  = [1e7, 1e3, 1e3]
γ0 = ωs ./ Qs

# Constants
ħ = 1.05457182e-34 * 1e9^2 / 1e6 # esta en nm i μs! 
kB = 1.380649e-23 * 1e9^2 / 1e6^2
T = 300#1e-3
const nth = 1/(exp(ħ*ωm/(kB*T))-1)


# Measurement "rates"
g0 = 2π*465 # Hz, pel mode principal. Busquem la g0 pels altres modes trobant la G: G = g0 / xzpf --> g0_altres = G * xzpf_altres
ncav = 1e5
xzp = @. sqrt(ħ/(2*m*ωs)) # en nm !
G = g0/xzp[1] # en Hz/nm
g0vec = G * xzp

g = g0vec*sqrt(ncav)
κ = 2π*45e6 # Hz
Γqba = 4*g.^2/κ /1e6 # en MHz?
η = 0.9

Γba_vec = Γqba
Γmeas_vec = Γqba*η


# Initial conditions and SDE parameters TODO: Change initial conditions to correct size
p = (γ0, ωm, ωs, Γba_vec ,Γmeas_vec)

VX_ss = @. -γ0 + sqrt(γ0^2 + 4*Γmeas_vec*(γ0)*(nth + 0.5) + Γba_vec) / (8*Γmeas_vec)
# In the case in which all 3 have same Q
C0 = [i == j ? VX_ss[1] : 0.0 for i in 1:6, j in 1:6]
u0 = vcat(zeros(6), vec(C0))

tspan = (0.0, 10000) # μs

prob= SDEProblem(moments_evolution_3modes, infogain_3modes, u0, tspan, p,  noise_rate_prototype = zeros(42, 2))
@time sol = solve(prob,
    SRA1(),
    saveat = 0.1/ωm,
    dt = 1e-3,
    dtmax = 1e-1,
    # dtmin = 1e-6,
    maxiters = tspan[2]*1e7,
    progress = true,
);

# Data visualization

pvars = plot(sol.t, sol[1,:], label = L"\langle X_1 \rangle", xlabel = "t [μs]", ylabel = "Nondimensional")
plot!(sol.t[:], sol[2,:], label = L"\langle Y_1 \rangle")
plot!(sol.t[:], sol[3,:], label = L"\langle X_2 \rangle")
plot!(sol.t[:], sol[4,:], label = L"\langle Y_2 \rangle")
plot!(sol.t[:], sol[5,:], label = L"\langle X_3 \rangle")
plot!(sol.t[:], sol[6,:], label = L"\langle Y_3 \rangle")

# Plotting time evolution of the variances (of the mode we aim to study, x1)
# pcovs = plot(sol.t, sol[7,:], label = L"V_{X_1}")
# plot!(sol.t, sol[8,:], label = L"C_{X_1Y_1}")
# plot!(sol.t, sol[9,:], label = L"C_{X_1X_2}")
# plot!(sol.t, sol[10,:], label = L"C_{X_1Y_2}")
# plot!(sol.t, sol[11,:], label = L"C_{X_1X_3}")
# plot!(sol.t, sol[12,:], label = L"C_{X_1Y_3}")

# plot!(sol.t, sol[14,:], label = L"V_{Y_1}", ls = :dash)
# plot!(sol.t, sol[15,:], label = L"C_{Y_1X_2}", ls = :dash)
# plot!(sol.t, sol[16,:], label = L"C_{Y_1Y_2}", ls = :dash)
# plot!(sol.t, sol[17,:], label = L"C_{Y_1X_3}", ls = :dash)
# plot!(sol.t, sol[18,:], label = L"C_{Y_1Y_3}", ls = :dash, xlabel = "t, μs", ylabel = "Covariances, adimensional")



# Looking at spectrum:
using DSP
# Aquests estan en hertz
κ1 = 10/11 * κ / 1e6
κ2 = 1/11 * κ / 1e6
Δ  = 0#-κ / (2*sqrt(3)) / 1e6

Q2 = [sol[1:2:5,i] .* cos.(ωs.*sol.t[i]) .+ sol[2:2:6,i] .* sin.(ωs.*sol.t[i]) for i in eachindex(sol.t)]

Q = [sol[2*i-1,:].* cos.(ωs[i].*sol.t) .+ sol[2i,:].* sin.(ωs[i].*sol.t) for i in 1:3]

Q1 = sol[1,:] .* cos.(ωs[1]*sol.t) .+ sol[2,:].*sin.(ωs[1]*sol.t)
Q2 = sol[3,:] .* cos.(ωs[2]*sol.t) .+ sol[4,:].*sin.(ωs[2]*sol.t)
Q3 = sol[5,:] .* cos.(ωs[3]*sol.t) .+ sol[6,:].*sin.(ωs[3]*sol.t)

δω = sum.([g0vec/1e6.*sol[1:2:5,i] for i in eachindex(sol.t)]) # aixo en MHz
Sin = 1

# psd = welch_pgram(sol[3,1:100:end]; onesided = false, fs = 0.01/ωm*100)
# plot!(psd.freq, psd.power, yscale = :log10)

S0 =  sqrt(κ1*κ2)/ (κ/(1e6)/2 + im*Δ)*Sin
S1 = -4im*sqrt(κ1*κ2)/ (κ/(1e6) + 2*im*Δ)^2*Sin
S2 = -16*sqrt(κ1*κ2)/ (κ/(1e6) + 2*im*Δ)^3*Sin

Sout = S0 .+ S1*δω .+ S2*(δω).^2
# Sout = S2*(δω).^2

i_signal = angle.(Sout)

# plot(i_signal)

# sol = sol[1:100:end]



using FFTW

fs = 1/(0.1/ωm)
F = fftshift(fft(Q[1]))
F2 = fftshift(fft(Q[2]))
F3 = fftshift(fft(Q[3]))
Fangle = fftshift(fft(i_signal))
freqs = fftshift(fftfreq(length(0:0.1/ωm:10000)+1, fs))

plot(freqs[1:5:end], abs.(Fangle[1:5:end]), yscale = :log10, xlims = (-1.5, 1.5))

plot(freqs, abs.(F), xlims = (-5,5), yscale = :log10)
plot!(freqs, abs.(F2), yscale = :log10)
plot!(freqs, abs.(F3), yscale = :log10)