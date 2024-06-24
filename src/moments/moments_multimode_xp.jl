# Script for time evolution of the moments and covariances (SDE system) including 3 modes.TODO: Add appropriate initial conditions
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
    # x, p, C
    xs = @view u[1:2:5]
    ps = @view u[2:2:6]
    C = @view u[7:end]
    # dx, dp, dC
    dxs = @view du[1:2:5]
    dps = @view du[2:2:6]
    dC = @view du[7:end]

    #First moments
    @. dxs = -γ0/2*xs + ωs*ps
    @. dps = -γ0/2*ps - ωs*xs

    # Covariance matrix
    for M = 1:2, N = 1:2
        for i in 1:3, j in 1:3
            # Getting the right indices
            Mi = right_index(M,i)
            Nj = right_index(N,j)
            Mconji = right_index(M,i,2)
            Nconjj = right_index(N,j,2)

            dC[colmaj(Mi, Nj)] = - (γ0[i]+γ0[j])/2 * C[colmaj(Mi,Nj)] + δ(Mi,Nj)*γ0[i]*(nth+0.5) + δ(M,2)*δ(N,2)*sqrt(Γba[i]*Γba[j]) +
                     + (-1)^δ(M,2)*ωs[i]*C[colmaj(Mconji, Nj)] + (-1)^δ(N,2)*ωs[j]*C[colmaj(Mi, Nconjj)] +
                     -4*(sum([sqrt(Γmeas[k])*C[colmaj(Mi, 2k-1)] for k in 1:3]))*(sum([sqrt(Γmeas[l])*C[colmaj(Nj, 2l-1)] for l in 1:3]))
        end
    end
    nothing
end

function infogain_3modes(du,u,p,t)
    @views Γmeas = p[end]
    C = @view u[7:end]
    du[1:2:5] = 2*[sum([sqrt(Γmeas[j])*C[colmaj(xi,2j-1)] for j in 1:3]) for xi in 1:2:5] # dWx in d<Xi>
    du[2:2:6] = 2*[sum([sqrt(Γmeas[j])*C[colmaj(pi,2j-1)] for j in 1:3]) for pi in 2:2:6] # dWx in d<Yi>
    du[7:end] .= 0.0
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
T = 300
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


# Initial conditions and SDE parameters TODO: Change initial conditions to correct steady state
p = (γ0, ωm, ωs, Γba_vec ,Γmeas_vec)

Vx_th = nth + 0.5
C0 = [i == j ? Vx_th : 0.0 for i in 1:6, j in 1:6]
u0 = vcat(zeros(6), vec(C0))

tspan = (0.0, 100) # μs

prob= SDEProblem(moments_evolution_3modes, infogain_3modes, u0, tspan, p)
@time sol = solve(prob,
    SOSRI(),
    saveat = 1/10,
    dt = 1e-14,
    dtmax = 1/10,
    dtmin = 1e-16,
    maxiters = tspan[2]*1e7,
    progress = true,
);

# Data visualization

pvars = plot(sol.t, sol[1,:], label = L"\langle x_1 \rangle", xlabel = "t [μs]", ylabel = L"X/x_{zpf}"*" (unitless)")
plot!(sol.t[:], sol[2,:], label = L"\langle p_1 \rangle", dpi = 600)
plot!(sol.t[:], sol[3,:], label = L"\langle x_2 \rangle")
plot!(sol.t[:], sol[4,:], label = L"\langle p_2 \rangle")
plot!(sol.t[:], sol[5,:], label = L"\langle x_3 \rangle")
plot!(sol.t[:], sol[6,:], label = L"\langle p_3 \rangle")

# Plotting time evolution of the variances (of the mode we aim to study, x1)
pcovs = plot(sol.t, sol[7,:], label = L"V_{x_1}")
plot!(sol.t, sol[8,:], label = L"C_{x_1p_1}")
plot!(sol.t, sol[9,:], label = L"C_{x_1x_2}")
plot!(sol.t, sol[10,:], label = L"C_{x_1p_2}")
plot!(sol.t, sol[11,:], label = L"C_{x_1x_3}")
plot!(sol.t, sol[12,:], label = L"C_{x_1p_3}")

plot!(sol.t, sol[14,:], label = L"V_{p_1}", ls = :dash)
plot!(sol.t, sol[15,:], label = L"C_{p_1x_2}", ls = :dash)
plot!(sol.t, sol[16,:], label = L"C_{p_1p_2}", ls = :dash)
plot!(sol.t, sol[17,:], label = L"C_{p_1x_3}", ls = :dash)
plot!(sol.t, sol[18,:], label = L"C_{p_1p_3}", ls = :dash, xlabel = "t, μs", ylabel = "Covariances, adimensional")



# Looking at spectrum:
# Aquests estan en hertz
κ1 = 10/11 * κ / 1e6
κ2 = 1/11 * κ / 1e6
Δ  = -κ / (2*sqrt(3)) / 1e6

δω = sum.([g0vec/1e6.*sol[1:2:5,i] for i in eachindex(sol.t)]) # aixo en MHz
Sin = 1

# psd = welch_pgram(sol[3,1:100:end]; onesided = false, fs = 0.01/ωm*100)
# plot!(psd.freq, psd.power, yscale = :log10)

S0 =  sqrt(κ1*κ2)/ (κ/(1e6)/2 + im*Δ)*Sin
S1 = -4im*sqrt(κ1*κ2)/ (κ/(1e6) + 2*im*Δ)^2*Sin
S2 = -16*sqrt(κ1*κ2)/ (κ/(1e6) + 2*im*Δ)^3*Sin

Sout = S0 .+ S1*δω .+ S2*(δω).^2

i_signal = angle.(Sout)



using FFTW

Q = sol

fs = 1/(0.1/ωm)
F = fftshift(fft(Q[1,:]))
F2 = fftshift(fft(Q[3,:]))
F3 = fftshift(fft(Q[5,:]))
Fangle = fftshift(fft(i_signal))
freqs = fftshift(fftfreq(length(0:1/10:tspan[2]), 1/10))


df = 10
plot_angle = plot(freqs[1:df:end], abs.(Fangle[1:df:end]), yscale = :log10, xlims = (-5, 5), xlabel = L"ω"*" [MHz]", label = L"arg(S_{out})", alpha = 1, dpi = 700, ylabel = "Fourier Transform, logarithmic (1/Hz)")

plot(freqs[1:df:end], abs.(F[1:df:end]), xlims = (-5,5), yscale = :log10, label = L"\mathcal{F}(x_1)", xlabel = "ω", alpha = 0.7)
plot!(freqs[1:df:end], abs.(F2[1:df:end]), yscale = :log10, label = L"\mathcal{F}(x_2)", alpha = 0.7)
plot!(freqs[1:df:end], abs.(F3[1:df:end]), yscale = :log10, label = L"\mathcal{F}(x_3)", alpha = 0.7) 
plot!(freqs[1:df:end], abs.(F[1:df:end]).+abs.(F2[1:df:end]).+abs.(F3[1:df:end]))