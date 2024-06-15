using FFTW

# Define the constants
const_val = 1.0  # You can change this value as needed

# Define the functions in the frequency domain
X(ω) = const_val * sinc(const_val * ω / (2 * pi))
H(ω) = exp.(-0.1 * abs(ω))

# Choose the number of samples (N) and the sampling rate (f_s)
N = 1024  # Number of samples, choose a power of two for efficiency
f_s = 1.0 # Sampling rate in Hz

# Define the frequency vector
frequencies = (0:N-1) * (f_s / N) .- (f_s / 2)  # Center the frequencies around zero
ω = 2 * pi * frequencies  # Convert frequencies to angular frequency

# Evaluate the functions at these frequency points
X_ω = X.(ω)
H_ω = H.(ω)

# Multiply the expressions in the frequency domain
Y_ω = X_ω .* H_ω

# Perform the inverse Fourier transform to get the time-domain signal
y_t = ifft(fftshift(Y_ω))  # Use fftshift to center the frequency domain data

# y_t now contains the time domain representation of the product of X(w) and H(w)
println(real(y_t))  # Print the real part of the result
plot(real(y_t))