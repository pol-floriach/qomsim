using DSP, Plots

fs = 1000

t = collect(0:1/fs:100)
x = sin.(2π*5*t)
y = sin.(2π*100*t)
z = x.+y

fc = 50/fs
order = 3
hp_filter = digitalfilter(Lowpass(fc), Butterworth(4))

filtered = filtfilt(hp_filter, z)

plot(t,x, xlims = (0,1))
plot!(t,y)
plot!(t,z)
plot!(t,filtered, lw = 2)
