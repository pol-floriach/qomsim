using PyCall, Plots, LaTeXStrings, FFTW, DSP
@pyimport scipy.signal as spsig
function my_psd(s,fs)
    N = length(s)
    Y = fft(s)
    P2 = abs.(Y / N)
    P1 = P2[1:N ÷ 2 + 1]
    P1[2:end-1] *= 2

    f = fs*(0:N ÷ 2) / N
    psd = P1.^2 / fs

    return f, psd
end
function theoretical_psd(Ω)
    χg = @. 1 / (1 + 2im*(Ω - Ω₀)/Γ₀)
    Sxx = abs2.(χg) * nth * 8 * xzpf^2 / Γ₀
    return Sxx
end

Ωvec = (0:.1:50) * 2π
Sxxt = theoretical_psd(Ωvec)


fs = 1/1e-2
f, Pxx_den = spsig.periodogram(sol_stoch[1, :], fs)
f2,psd = my_psd(sol_stoch[1,:],fs)
psd2 = welch_pgram(sol_stoch[1,:]; fs = fs)


p2 = plot(f[1:10:end], Pxx_den[1:10:end],
    yscale=:log10, 
    xlabel="frequency = [MHz]", 
    ylabel="PSD [1/MHz] (log)", 
    label = "Scipy periodogram",
    title=L"S_{xx}",
    markersize = 1,
    markerstrokewidth = 0,
    ylims = (1e-20, 1),
    
)
plot!(psd2.freq, psd2.power, label = "DSP.jl welch")
plot(Ωvec/2π, Sxxt, yscale = :log10, lw = 2, label = "Theoretical")