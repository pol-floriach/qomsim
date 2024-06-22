# Script for time evolution of the moments and covariances (SDE system) including 3 modes and filtered feedback.
using DifferentialEquations, Plots, LaTeXStrings, ProgressLogging, DelimitedFiles

# Function to obtain O_k index of matrix (in a row or column)
function right_index(O,k,target= 1)
    return O == target ?  2k-1 : 2k
end
# Function to get index of a matrix stored in an array in column-major order (for 6x6 matrix)
colmaj(row, col) = @. row + (col-1)*4

# Kronecker delta function for better legibility
δ(x,y) = ==(x,y)

function shift_push!(meas_xs, new)
    popfirst!(meas_xs)
    push!(meas_xs, new)
end

function apply_filter(meas_xs, coeffs)
    sum(meas_xs.*coeffs)
end

function apply_iir_filter(meas_xs, meas_filtered,coeffs_num, coeffs_den)
    1/coeffs_den[1] * sum(coeffs_num.*meas_xs) - sum(coeffs_den.*meas_filtered)
end

# Time evolution w/o information gain. TODO If too many allocations: try with staticarray / add Mi, Nj etc as parameters and mutate
function moments_evolution_2modes(du, u, p, t)
    # Parameters
    @views γ0, ωm, ωs, Γba ,Γmeas, tsamp, meas_buffer, _, _, _, _, _, _ = p
    # x, p, C
    xs = @view u[1:2:3]
    ps = @view u[2:2:4]
    C = @view u[5:end]
    # dx, dp, dC
    dxs = @view du[1:2:3]
    dps = @view du[2:2:4]
    dC = @view du[5:end]

    #First moments
    @. dxs = -γ0/2*xs + ωs*ps
    @. dps = -γ0/2*ps - ωs*xs -  F_fbvec_low[end] - F_fbvec_high[end]
    
    # Covariance matrix
    for M = 1:2, N = 1:2
        for i in 1:2, j in 1:2
            # Getting the right indices
            Mi = right_index(M,i)
            Nj = right_index(N,j)
            Mconji = right_index(M,i,2)
            Nconjj = right_index(N,j,2)

            dC[colmaj(Mi, Nj)] = - (γ0[i]+γ0[j])/2 * C[colmaj(Mi,Nj)] + δ(Mi,Nj)*γ0[i]*(nth+0.5) + δ(M,2)*δ(N,2)*sqrt(Γba[i]*Γba[j]) +
                     + (-1)^δ(M,2)*ωs[i]*C[colmaj(Mconji, Nj)] + (-1)^δ(N,2)*ωs[j]*C[colmaj(Mi, Nconjj)] +
                     -4*(sum([sqrt(Γmeas[k])*C[colmaj(Mi, 2k-1)] for k in 1:2]))*(sum([sqrt(Γmeas[l])*C[colmaj(Nj, 2l-1)] for l in 1:2]))
        end
    end
    nothing
end

function infogain_2modes(du,u,p,t)
    @views Γmeas = p[5]
    C = @view u[7:end]
    du[1:2:3] = 2*[sum([sqrt(Γmeas[j])*C[colmaj(xi,2j-1)] for j in 1:2]) for xi in 1:2:3] # dWx in d<Xi>
    du[2:2:4] = 2*[sum([sqrt(Γmeas[j])*C[colmaj(pi,2j-1)] for j in 1:2]) for pi in 2:2:4] # dWx in d<Yi>
    du[5:end] .= 0.0
end

function affect!(int)
    shift_push!(int.p[7], sum(int.u[1:2:3]))
    F_fb_low =      apply_filter(int.p[7][int.p[10]], int.p[8])
    F_fb_high = 1.5*apply_filter(int.p[7][int.p[11]], int.p[9])
    push!(int.p[12], F_fb_low)
    push!(int.p[13], F_fb_high)
end

# System parameters
begin
    # Parameters
    m = 1e-12
    ωm = 2π*1.1
    ωs = [2π*0.5, 2π*1.5]
    Qs  = [1e3, 1e3]
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
end

    # Parameters for the filter
begin
    fsamp = 125 # Sampling frequency, MHz
    tsamp = 1/fsamp # μs

    period_low = 2π / ωs[1] 
    period_4_low = period_low/4

    period_high = 2π / ωs[2]
    period_4_high = period_high/4

    # period_band = 2π / ωs[1]
    # period_4_band = period_band/4

# Filter parameters (for low and high pass)

    N_filter = 30+1
    # N_band =   85+1

    Nfix_low =  ceil(Int, 0.85 / tsamp)
    Nfix_high = ceil(Int, 0.21 / tsamp)
    # Nfix_band = ceil(Int, 0.1  / tsamp)

    NT4_low = ceil(Int,  period_4_low  / tsamp)
    NT4_high = ceil(Int, period_4_high / tsamp)
    # NT4_band = ceil(Int, period_4_band / tsamp)

    N_buffer_low = (N_filter  +  NT4_low  + Nfix_low)
    N_buffer_high = (N_filter +  NT4_high + Nfix_high)
    # N_buffer_band = (N_band   + NT4_band  + Nfix_band)
    N_buffer = max(N_buffer_high, N_buffer_low)
    meas_buffer = zeros(N_buffer)
    F_fb = 0.0

    range_low = 1:N_filter
    range_high = 1+Nfix_high: N_filter + Nfix_high
    # range_band = 1+Nfix_band: N_band + Nfix_band

    # Reading .dat files containing filter coefficients
    coeffs_lp = vec(readdlm("lowpass_rectangular.dat"))
    coeffs_hp = vec(readdlm("highpass_rectangular.dat"))
    # coeffs_bp = vec(readdlm("bandpass_rectangular.dat"))
    F_fbvec_low = [F_fb]
    F_fbvec_high = [F_fb]
    # F_fbvec_band = [F_fb]
end

# Initial conditions and SDE parameters 
# p = (γ0, ωm, ωs, Γba_vec ,Γmeas_vec, tsamp, meas_buffer, coeffs_lp, coeffs_hp, coeffs_bp, range_low, range_high, range_band, F_fbvec_low, F_fbvec_high, F_fbvec_band)
p = (γ0, ωm, ωs, Γba_vec ,Γmeas_vec, tsamp, meas_buffer, coeffs_lp, coeffs_hp, range_low, range_high, F_fbvec_low, F_fbvec_high)

begin
    Vx_th = nth + 0.5;
    C0 = [i == j ? Vx_th : 0.0 for i in 1:4, j in 1:4];
    u0 = vcat(zeros(4), vec(C0));
end
tspan = (0.0, 50) # μs

measuring_times = collect(0:tsamp:tspan[2])
measuring_times = measuring_times[3000:end]
cb = PresetTimeCallback(measuring_times, affect!)

prob= SDEProblem(moments_evolution_2modes, infogain_2modes, u0, tspan, p)
@time sol = solve(prob,
    SOSRI(),
    saveat = tsamp,
    dt = 1e-14,
    dtmax = tsamp,
    dtmin = 1e-16,
    maxiters = tspan[2]*1e7,
    progress = true,
    save_idxs = 1:6,
    callback = cb
);

using FFTW
Fforce = fftshift(fft(F_fbvec_band[2:end]))
freqsforce = fftshift(fftfreq(length(measuring_times), fsamp))
plot!(freqsforce, abs.(Fforce), xlims = (-2.5,2.5))


# Data visualization

pvars = plot(sol.t, sol[1,:], label = L"\langle x_1 \rangle", xlabel = "t [μs]", ylabel = "Nondimensional")
# plot!(sol.t[:], sol[2,:], label = L"\langle p_1 \rangle")
plot!(sol.t[:], sol[3,:], label = L"\langle x_2 \rangle")
# plot!(sol.t[:], sol[4,:], label = L"\langle p_2 \rangle")


plot!(measuring_times, F_fbvec_low[2:end], label = "Filter lp", xlabel = "t [μs]", ylabel = "Nondimensional")
plot!(measuring_times, F_fbvec_high[2:end], label = "Filter hp", xlabel = "t [μs]", ylabel = "Nondimensional")
plot!(measuring_times, F_fbvec_band[2:end], label = "Filter bp", xlabel = "t [μs]", ylabel = "Nondimensional")

# Plotting time evolution of the variances (of the mode we aim to study, x1)
# pcovs = plot(sol.t, sol[7,:], label = L"V_{x_1}")
# plot!(sol.t, sol[8,:], label = L"C_{x_1p_1}")
# plot!(sol.t, sol[9,:], label = L"C_{x_1x_2}")
# plot!(sol.t, sol[10,:], label = L"C_{x_1p_2}")
# plot!(sol.t, sol[11,:], label = L"C_{x_1x_3}")
# plot!(sol.t, sol[12,:], label = L"C_{x_1p_3}")

# plot!(sol.t, sol[14,:], label = L"V_{p_1}", ls = :dash)
# plot!(sol.t, sol[15,:], label = L"C_{p_1x_2}", ls = :dash)
# plot!(sol.t, sol[16,:], label = L"C_{p_1p_2}", ls = :dash)
# plot!(sol.t, sol[17,:], label = L"C_{p_1x_3}", ls = :dash)
# plot!(sol.t, sol[18,:], label = L"C_{p_1p_3}", ls = :dash, xlabel = "t, μs", ylabel = "Covariances, adimensional")



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

# fs = 1/(0.1/ωm)
F = fftshift(fft(Q[1,:]))
F2 = fftshift(fft(Q[3,:]))
F3 = fftshift(fft(Q[5,:]))

F =  fftshift(fft(Q[1,Int(ceil(10*end/11+1)):end]))
F2 = fftshift(fft(Q[3,Int(ceil(10*end/11+1)):end]))
F3 = fftshift(fft(Q[5,Int(ceil(10*end/11+1)):end]))
Fangle = fftshift(fft(i_signal))
freqs = fftshift(fftfreq(length(0:tsamp:tspan[2]), fsamp))
freqs2 = fftshift(fftfreq(length(F3), fsamp))


df = 1
# plot_angle = plot!(freqs[1:df:end], abs.(Fangle[1:df:end]), yscale = :log10, xlims = (-5, 5), xlabel = L"ω"*" [MHz]", label = L"arg(S_{out})", alpha = 1, dpi = 700, ylabel = "Fourier Transform, logarithmic (1/Hz)")

plot(freqs2[1:df:end], abs.(F[1:df:end]), xlims = (-5,5), yscale = :log10, label = L"\mathcal{F}(x_1), fb", xlabel = L"\omega"*" [MHz]", alpha = 0.7, ylabel = "Spectra [1/√MHz]")
plot!(freqs2[1:df:end], abs.(F2[1:df:end]), yscale = :log10, label = L"\mathcal{F}(x_2), fb", alpha = 0.7)
plot!(freqs2[1:df:end], abs.(F3[1:df:end]), yscale = :log10, label = L"\mathcal{F}(x_3), fb", alpha = 0.7) 