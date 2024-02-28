# Note: Work in progress
using Distributions
function kalman_step(zk, xpk, Pk, Φ, Q, H, R)
    Kk = Pk * transpose(H) * inv(H*Pk*transpose(H) + R)
    xk = xpk + Kk * (zk - H * xpk)
    Pk = (I - Kk * H) * Pk
    xpk1 = Φ * Pk * transpose(Phi) + Q

    return xk, xpk1, Pkp1, Kk
end

H = [1 0; 0 0]
R = [1 0; 0 1e15]
xps = zeros(length(t),2) # Position estimates

Ps = zeros(length(t),size(R)...)
Ps[1,:,:] = R

# Simulating measurement values (change for actual measurements)
zs = [zeros(2) for i in 1:length(t)]
zs[1] = H * xs[1] + rand(MvNormal(zeros(2),R))

xps[1][1,1] = zs[1][1,1]
xps[1][2] = 0

xfs = [zeros(2) for i in 0:length(t)-1]

Ks = zeros(length(t),size(R)...)