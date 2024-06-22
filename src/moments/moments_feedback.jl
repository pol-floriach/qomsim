# Script for time evolution of the moments and variances (SDE system) with feedback
using DifferentialEquations, StaticArrays, Plots, LaTeXStrings, ProgressLogging, DelimitedFiles

# Moments time evolution (deterministic part) with feedback term in Hamiltonian: H_fb(t) = -G⋅⟨x(t-T/4)⟩x(t)
function moments_evolution(u,p,t)
    @views γ0, ω, μ, η, _, _, _, _, _, _, _, _, F_fbvec, _, _, _, _, _, _, _, _, _ = p

    @views x, p, Vx, Vp, Cxp = u

    dx    = -γ0/2*x + ω*p
    dp    = -γ0/2*p - ω*x  - F_fbvec[end]


    dVx   = -γ0*Vx + 2*ω*Cxp - 4*μ*η*Vx^2  + γ0*(nth+0.5)
    dVp   = -γ0*Vp - 2*ω*Cxp - 4*μ*η*Cxp^2 + γ0*(nth+0.5) + μ
    dCxp  = ω*Vp - ω*Vx - γ0*Cxp - 4*μ*η*Vx*Cxp
    SA[dx, dp, dVx, dVp, dCxp]
end

# Stochastic part of the equations (due to information gain)
function moments_infogain(u,p,t)
     μ, η = p[3:4]
     Vx, Cxp = @view u[3:2:5]
    sqrt(4*μ*η)*SA[Vx, Cxp, 0, 0, 0]
end

function affect!(int)
    # F_fb = apply_filter(int.p[6][1:int.p[7]], int.p[8])
    _, _, _, _, _, meas_buffer, Nfilterlp, Nfilterbp, Nfilterhp, coeffs_lp, coeffs_hp, coeffs_bp, F_fbvec, Ndelay_lp, Ndelay_bp, Ndelay_hp, u_buffer_lp, u_buffer_bp, u_buffer_hp, modulated_meas_buffer, u_buffer_bp_mod, F_fbvec_mod = int.p

    shift_push!(meas_buffer, int.u[1])
    # shift_push!(modulated_meas_buffer, int.u[1]*sin(ωmid*int.t))

    F_fb_lp = apply_iir_filter(meas_buffer[1+Ndelay_lp:Ndelay_lp+Nfilterlp], u_buffer_lp[1:Nfilterlp-1],coeffs_lp[1:Nfilterlp], coeffs_lp[Nfilterlp+1:end])
    F_fb_bp = apply_iir_filter(meas_buffer[1+Ndelay_bp:Ndelay_bp+Nfilterbp], u_buffer_bp[1:Nfilterbp-1],coeffs_bp[1:Nfilterbp], coeffs_bp[Nfilterbp+1:end])
    F_fb_hp = apply_iir_filter(meas_buffer[1+Ndelay_hp:Ndelay_hp+Nfilterhp], u_buffer_hp[1:Nfilterhp-1], coeffs_hp[1:Nfilterhp], coeffs_hp[Nfilterhp+1:end])
    
    # F_fb_bp = bandpass_demodulating(modulated_meas_buffer[1+Ndelay_bp:Ndelay_bp+Nfilter1], u_buffer_bp_mod[1:Nfilter1-1], ωmid*(int.t-tau), coeffs_bp[1:Nfilter1], coeffs_bp[Nfilter1+1:end], F_fbvec_mod)
    push!(F_fbvec, 0.0*F_fb_lp + 1*F_fb_bp + 0.0*F_fb_hp)

    shift_push!(u_buffer_lp, F_fb_lp)
    shift_push!(u_buffer_bp, F_fb_bp)
    shift_push!(u_buffer_hp, F_fb_hp)
end

function shift_push!(meas_xs, new)
    pop!(meas_xs)
    pushfirst!(meas_xs, new)
end

function apply_fir_filter(meas_xs, coeffs)
    sum(meas_xs.*coeffs)
end

function apply_iir_filter(meas_xs, meas_filtered,coeffs_num, coeffs_den)
    sum_zeros = sum(coeffs_num.*meas_xs)
    sum_poles = sum(coeffs_den[2:end].*meas_filtered)
    1/coeffs_den[1] * (sum_zeros - sum_poles)
end

function bandpass_demodulating(signal_modulated, f_fb_modulatedbuffer, ωt, coeffs_num, coeffs_den, F_fbvec_mod)
    F_fb_modulated = apply_iir_filter(signal_modulated, f_fb_modulatedbuffer, coeffs_num, coeffs_den)
    shift_push!(f_fb_modulatedbuffer, F_fb_modulated)
    push!(F_fbvec_mod, F_fb_modulated)
    return F_fb_modulated #* sin(ωt)
end

# Parameters
begin
    ω = 2π*1.5# [ω] = MHz
    const ωmid = ω
    Q = 1e8
    γ0 = ω / Q

    # Measurement "rate" ([MHz])
    g0 = 2π*465
    ncav = 1e5

    g = g0*sqrt(ncav)
    κ = 2π*45e6
    μ = 4*g^2/κ /1e6
    η = 0.9


    # Constants
    const ħ = 1.05457182e-34 * 1e9^2 / 1e6
    const kB = 1.380649e-23 * 1e9^2 / 1e6^2
    T = 300
    const nth = 1/(exp(ħ*ω/(kB*T))-1)
    Vx_th = Vp_th = nth + 0.5
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
    # coeffs_hp = vec(readdlm("filter_coeffs/highpass_iir_6.dat"))
    coeffs_hp = vec(readdlm("filter_coeffs/bandpass_iir_10_alt.dat"))
    # coeffs_bp = vec(readdlm("filter_coeffs/lowpass_iir_4_pelbandpass.dat"))

    F_fbvec = Float64[]
    F_fbvec_mod = Float64[]
    ubuffer_hp = zeros(N_buffer-1)
    ubuffer_bp = zeros(N_buffer-1)
    ubuffer_lp = zeros(N_buffer-1)

    modulated_meas_buffer = zeros(N_buffer)
    u_buffer_bp_mod = zeros(N_buffer-1)
end



tau = 0.0
# Initial conditions and SDE parameters
p = (γ0, ω, μ, η, tsamp, meas_buffer, Nfilterlp, Nfilterbp, Nfilterhp, coeffs_lp, coeffs_hp, coeffs_bp, F_fbvec, Ndelay_lp, Ndelay_bp, Ndelay_hp, ubuffer_lp, ubuffer_bp, ubuffer_hp, modulated_meas_buffer, u_buffer_bp_mod, F_fbvec_mod,)
u0 = [0, 0, Vx_th, Vp_th, 0.0]
tspan = (0.0, 100) # μs

measuring_times = collect(0:tsamp:tspan[2])
cb = PresetTimeCallback(measuring_times, affect!, save_positions = (false, false))
# Simulation
prob= SDEProblem(moments_evolution, moments_infogain, u0, tspan, p)
@time sol = solve(prob,
    SOSRI(),
    saveat = tsamp,
    dtmin = 1e-16,
    dtmax = tsamp,
    dt = 1e-15,
    maxiters = tspan[2]*1e7,
    progress = true,
    callback = cb
);

pmeans = plot(sol.t, sol[1,:], label = L"\langle x \rangle", xlabel = "t [μs]", ylabel = "Nondimensional")
plot!(sol.t, sol[2,:], label = L"\langle p \rangle")
plot!(measuring_times, F_fbvec , label = "Filter", xlabel = "t [μs]", ylabel = "Nondimensional",xlims = (10,13))



function energy(sol)
    return 1/2 * (sol[1,:].^2 + sol[2,:].^2)
end

Es = energy(sol)
Es_ratio1 = Es / (3/2)

es = plot!(sol.t[2:end], Es_ratio1[1,2:end], label = L"E_1")
plot!(sol.t[2:end], Es_ratio1[2,2:end], label = L"E_2")
plot!(sol.t[2:end], Es_ratio1[3,2:end], label = L"E_3", yscale = :log10, xlabel = "t [μs]", ylabel = "Energy ratio", dpi = 500)







# using FFTW
# Fforce = fftshift(fft(F_fbvec[2:end]))
# freqsforce = fftshift(fftfreq(length(measuring_times), fsamp))
# plot(freqsforce, abs.(Fforce), xlims = (-2.5,2.5))



# pvars = scatter(sol.t, sol[,:], label = L"V_x", xlabel = "t [μs]", ylabel = "Nondimensional")
# plot!(sol.t, sol[3,:], label = L"V_p")
