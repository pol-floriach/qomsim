using DSP, Plots

fs = 125

ω1 = 2π * 0.5
ω2 = 2π * 1.6

t = collect(0:1/fs:100)
x = sin.(ω1*t)
y = sin.(ω2*t)
z =  x .+y

fc = (1.1 + 0.25)/(fs/2)
fc2 = (1.1 - 0.25)/(fs/2)
fc3 = (1.1 - 0.01)/(fs/2)
fc4 = (1.1 + 0.01)/(fs/2)



hp_filter = digitalfilter(Highpass(fc), Butterworth(10))
lp_filter = digitalfilter(Lowpass(fc2), Butterworth(10))
bp_filter = digitalfilter(Bandpass(fc3,fc4), Butterworth(10))
# hp_filter = digitalfilter(Highpass(fc;fs = fs), FIRWindow(rect(21)))


filtered = filt(hp_filter, x)
filtered2 = filt(lp_filter, x)


plot(t,x, xlims = (0,10))
plot(t,y)
plot!(t .+ .0,filtered, lw = 2, xlims = (5,10))
plot!(t .+ .0,filtered2, lw = 2, xlims = (90,100))


hp_polynomial = (convert(PolynomialRatio, hp_filter))
num_coefs = coefb(hp_polynomial)
den_coefs = coefa(hp_polynomial)

lp_polynomial = (convert(PolynomialRatio, hp_filter))
num_coefs = coefb(hp_polynomial)
den_coefs = coefa(hp_polynomial)

bp_polynomial = (convert(PolynomialRatio, bp_filter))
num_coefs = coefb(bp_polynomial)
den_coefs = coefa(bp_polynomial)


# Provant la meva funcio (amb lowpass ara)
begin
    fsamp = 125 # Sampling frequency, MHz
    tsamp = 1/fsamp # μs

    period = 2π / ω1
    period_4 = period/4
end
begin
    N_filter = 10+1
    Nfix = ceil(Int, 0 / tsamp)
    # N_buffer = (N_filter +  ceil(Int, period_4 / tsamp)+ Nfix)
    N_buffer = N_filter

    meas_buffer = zeros(N_buffer)

    F_fbvec = Float64[]
    u_buffer = zeros(N_buffer-1)
end

function shift_push!(vector, new)
    popfirst!(vector)
    push!(vector, new)
end

function apply_iir_filter(meas_xs, meas_filtered,num_coefs, den_coefs)
    1/den_coefs[1] * (sum(reverse(num_coefs).*meas_xs) - sum(reverse((den_coefs[2:end])).*meas_filtered))
end

function affect!(u,p)
    shift_push!(p[6], u[1])
    # if i >= N_filter
    F_fb = apply_iir_filter(p[6][1:p[7]], p[11][1:p[7]-1], p[8][1:11], p[8][12:22])
    push!(p[9], F_fb)
    shift_push!(p[11], F_fb)
    # end
end

# for i in eachindex(t)
#     # Add measurement x[n] to the buffer
#     shift_push!(meas_buffer, x[i])
#     # Calculate y[n]
#     if i >= N_filter
#     F_fb = apply_iir_filter(meas_buffer[1:N_filter], u_buffer[1:N_filter-1], num_coefs,den_coefs)
#     # Add it to the filter output record
#     push!(F_fbvec, F_fb)
#     # Add it to the output buffer vector
#     shift_push!(u_buffer, F_fb)
#     end
# end

coeffs_hp = vcat(num_coefs, den_coefs)

p = (γ0, ω, μ, η, tsamp, meas_buffer, N_filter, coeffs_hp, F_fbvec, Nfix, u_buffer)

for i in eachindex(t)
    affect!(x[i],p) 
end


# plot(t,y, xlims = (0,1))
plot(t, F_fbvec)

sum(F_fbvec.^2) / length(F_fbvec)