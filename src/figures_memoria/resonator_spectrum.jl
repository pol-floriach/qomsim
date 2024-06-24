using DelimitedFiles, Plots

data= readdlm("PSD_ample.csv", ',')[2:end,:]

freq = vec(data[:,1])./1e6
PSD = vec(data[:,2])

plot(freq, PSD, 
    legend = false, 
    xlabel = "f [MHz]", 
    ylabel = "PSD, arbitrary units",
    dpi = 600,
    xlims = (0,2),
    ylims = (-100,-40))

savefig("psd_resonator.svg")