# Script for time evolution of the moments and covariances (SDE system) including 3 modes and filtered feedback.
using DifferentialEquations, Plots, LaTeXStrings, ProgressLogging, DelimitedFiles

# Function to obtain O_k index of matrix (in a row or column)
function right_index(O,k,target= 1)
    return O == target ?  2k-1 : 2k
end
# Function to get index of a matrix stored in an array in column-major order (for 6x6 matrix)
colmaj(row, col) = @. row + (col-1)*6

# Kronecker delta function for better legibility
δ(x,y) = ==(x,y)

function shift_push!(meas_xs, new)
    pop!(meas_xs)
    pushfirst!(meas_xs, new)
end


function apply_iir_filter(meas_xs, meas_filtered,coeffs_num, coeffs_den)
    sum_zeros = sum(coeffs_num.*meas_xs)
    sum_poles = sum(coeffs_den[2:end].*meas_filtered)
    1/coeffs_den[1] * (sum_zeros - sum_poles)
end


# Time evolution w/o information gain. TODO If too many allocations: try with staticarray / add Mi, Nj etc as parameters and mutate
function moments_evolution_3modes(du, u, p, t)
    # Parameters
    @views γ0, ω, ωs, Γba, Γmeas, tsamp, _, _, _, _, _, _, _, _, _, _, _, _, _, _, _, _, _ = p
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
    @. dps = -γ0/2*ps - ωs*xs - F_fbvec[end] 
    
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
    @views Γmeas = p[5]
    C = @view u[7:end]
    du[1:2:5] = 2*[sum([sqrt(Γmeas[j])*C[colmaj(xi,2j-1)] for j in 1:3]) for xi in 1:2:5] # dWx in d<Xi>
    du[2:2:6] = 2*[sum([sqrt(Γmeas[j])*C[colmaj(pi,2j-1)] for j in 1:3]) for pi in 2:2:6] # dWx in d<Yi>
    du[7:end] .= 0.0
end


function affect!(int)
    _, _, _, _, _, _, meas_buffer, Nfilterlp, Nfilterbp, Nfilterhp, coeffs_lp, coeffs_hp, coeffs_bp, F_fbvec, Ndelay_lp, Ndelay_bp, Ndelay_hp, u_buffer_lp, u_buffer_bp, u_buffer_hp, modulated_meas_buffer, u_buffer_bp_mod, F_fbvec_mod = int.p

    shift_push!(meas_buffer, sum(int.u[1:2:5]))

    F_fb_lp = apply_iir_filter(meas_buffer[1+Ndelay_lp:Ndelay_lp+Nfilterlp], u_buffer_lp[1:Nfilterlp-1],coeffs_lp[1:Nfilterlp], coeffs_lp[Nfilterlp+1:end])
    F_fb_bp = apply_iir_filter(meas_buffer[1+Ndelay_bp:Ndelay_bp+Nfilterbp], u_buffer_bp[1:Nfilterbp-1],coeffs_bp[1:Nfilterbp], coeffs_bp[Nfilterbp+1:end])
    F_fb_hp = apply_iir_filter(meas_buffer[1+Ndelay_hp:Ndelay_hp+Nfilterhp], u_buffer_hp[1:Nfilterhp-1], coeffs_hp[1:Nfilterhp], coeffs_hp[Nfilterhp+1:end])
    
    # F_fb_bp = bandpass_demodulating(modulated_meas_buffer[1+Ndelay_bp:Ndelay_bp+Nfilter1], u_buffer_bp_mod[1:Nfilter1-1], ωmid*(int.t-tau), coeffs_bp[1:Nfilter1], coeffs_bp[Nfilter1+1:end], F_fbvec_mod)
    push!(F_fbvec, 0.5*F_fb_lp + 0.4*F_fb_bp + 0.5*F_fb_hp)

    shift_push!(u_buffer_lp, F_fb_lp)
    shift_push!(u_buffer_bp, F_fb_bp)
    shift_push!(u_buffer_hp, F_fb_hp)
end

function energy(sol)
    return 1/2 * (sol[1:2:5,:].^2 + sol[2:2:6,:].^2)
end




# System parameters
begin
    # Parameters
    m = 1e-12
    ωm = 2π*1.1
    ωs = [2π*1.1, 2π*1.5, 2π*0.5]
    Qs  = [1e8, 1e3, 1e3]
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

begin
    fsamp = 125 # Sampling frequency, MHz
    tsamp = 1/fsamp # μs

    period = 1 / 0.5
    periodlp_4 = period/4

    period = 1 / 1.1
    periodbp_4 = period/4

    period = 1 / 1.5
    periodhp_4 = period/4
end

# Filter parameters
begin
    Nfilterlp = 6+1
    Nfilterbp = 8+1
    Nfilterhp = 10+1

    Ndelay_hp = ceil(Int, fsamp*(0.11))
    Ndelay_bp = ceil(Int, fsamp*(0.66))
    Ndelay_lp = ceil(Int, fsamp*(periodlp_4 + 0.48))

    N_buffer = 200
    meas_buffer = zeros(N_buffer)
    # Reading .dat files containing filter coefficients
    coeffs_lp = vec(readdlm("filter_coeffs/lowpass_iir_6.dat"))
    coeffs_bp = vec(readdlm("filter_coeffs/bandpass_iir_8.dat"))
    coeffs_hp = vec(readdlm("filter_coeffs/bandpass_iir_10_alt.dat"))

    F_fbvec = Float64[]
    F_fbvec_mod = Float64[]
    ubuffer_hp = zeros(N_buffer-1)
    ubuffer_bp = zeros(N_buffer-1)
    ubuffer_lp = zeros(N_buffer-1)

    modulated_meas_buffer = zeros(N_buffer)
    u_buffer_bp_mod = zeros(N_buffer-1)
end

# Initial conditions and SDE parameters 
p = (γ0, ωm, ωs, Γba_vec, Γmeas_vec, tsamp, meas_buffer, Nfilterlp, Nfilterbp, Nfilterhp, coeffs_lp, coeffs_hp, coeffs_bp, F_fbvec, Ndelay_lp, Ndelay_bp, Ndelay_hp, ubuffer_lp, ubuffer_bp, ubuffer_hp, modulated_meas_buffer, u_buffer_bp_mod, F_fbvec_mod,)

begin
    Vx_th = nth + 0.5;
    C0 = [i == j ? Vx_th : 0.0 for i in 1:6, j in 1:6];
    u0 = vcat(zeros(6), vec(C0));
end

tspan = (0.0, 300) # μs

measuring_times = collect(0:tsamp:tspan[2])
cb = PresetTimeCallback(measuring_times, affect!)

prob= SDEProblem(moments_evolution_3modes, infogain_3modes, u0, tspan, p)
@time sol = solve(prob,
    saveat = tsamp,
    dt = 1e-14,
    dtmax = tsamp,
    dtmin = 1e-16,
    maxiters = tspan[2]*1e7,
    progress = true,
    save_idxs = 1:6,
    callback = cb
);

using DataFrames, CSV
df = DataFrame(sol)
rename!(df, [:timestep, :xmain, :pmain, :xlow, :plow, :xhigh, :phigh])

CSV.write("3modes_gain_04.csv", df)

# Data visualization

pvars = plot(sol.t, sol[1,:], label = L"\langle x_1 \rangle", xlabel = "t [μs]", ylabel = "Nondimensional")
# plot!(sol.t[:], sol[2,:], label = L"\langle p_1 \rangle")
plot!(sol.t[:], sol[3,:], label = L"\langle x_2 \rangle")
# plot!(sol.t[:], sol[4,:], label = L"\langle p_2 \rangle")
plot!(sol.t[:], sol[5,:], label = L"\langle x_3 \rangle")
# plot!(sol.t[:], sol[6,:], label = L"\langle p_3 \rangle")

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



Es = energy(sol)
Es_ratio1 = Es / (3/2)

stopplot = 35000
es = plot(sol.t[2:stopplot], Es_ratio1[1,2:stopplot], label = L"E_{main}")
plot!(sol.t[2:stopplot], Es_ratio1[2,2:stopplot], label = L"E_{low}")
plot!(sol.t[2:stopplot], Es_ratio1[3,2:stopplot], label = L"E_{high}", yscale = :log10, xlabel = "t [μs]", ylabel = "Energy ratio (unitless)", dpi = 500, size = (600, 350))


using LsqFit
@. model(x,p) = p[1] * x + p[2]

fit13 = curve_fit(model, sol.t[2:15000], log10.(Es_ratio1[1,2:15000]), [-10, 250.0])

es = plot(sol.t[2:stopplot], log10.(Es_ratio1[1,2:stopplot]), label = L"E_{main}")
plot!(sol.t[2:stopplot], log10.(Es_ratio1[2,2:stopplot]), label = L"E_{low}")
plot!(sol.t[2:stopplot], log10.(Es_ratio1[3,2:stopplot]), label = L"E_{high}", xlabel = "t [μs]", ylabel = "Energy ratio (unitless, logarithmic)", dpi = 500, size = (600, 350))

plot!(sol.t[2:stopplot], coef(fit13)[1].*sol.t[2:stopplot] .+ coef(fit13)[2], label = "")



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
