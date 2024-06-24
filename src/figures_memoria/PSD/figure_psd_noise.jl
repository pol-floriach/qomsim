using DelimitedFiles, Plots, LaTeXStrings

freqs = vec(readdlm("freqs_Hz"))
a = vec(readdlm("goodunits_PSD_Teff16mK"))
b = vec(readdlm("goodunits_PSD_Teff16mK_fit"))
c = vec(readdlm("goodunits_Sxx_Teff16mK_fit"))

plot(freqs,a, yscale = :log10, xlabel = "Frequency [Hz]", ylabel = "PSD [m²/Hz]", label = "Experiment")
plot!(freqs,b, label = L"S_{XX} + S_{imp}")
plot!(freqs, c, label = L"S_{XX}", dpi = 600)

