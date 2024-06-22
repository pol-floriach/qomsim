using DSP, Plots

fs = 2000

t = collect(0:1/fs:100)
x = sin.(2π*5*t)
y = sin.(2π*100*t)
z = x.+y
z2 = sin.(2π*5*t .- π/10)

fc = 6/fs


hp_filter = digitalfilter(Highpass(fc), Butterworth(6))
hp_filter = digitalfilter(Highpass(fc;fs = fs), FIRWindow(rect(21)))


filtered = filtfilt(hp_filter, z)

plot(t,x, xlims = (0,1))
plot(t,y);
plot!(t,z)
plot!(t,filtered, lw = 2, xlims = (0,1))

h = remez(151, )