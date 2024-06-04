using DifferentialEquations, Plots, LaTeXStrings

# Function to get matrix stored in an array in column-major order
getindex(row, col) = row + (col-1)*6
# Kronecker delta function for better legibility
δ(x,y) = ==(x,y)
# Time evolution w/o information gain TODO: canviar parametres als que toca
function moments_evolution_3modes(du, u, p, t)
    # Parameters
    @views γ0, ω, k, η = p
    # X, Y, C
    Xs = @view u[1:2:5]
    Ys = @view u[2:2:6]
    C = @view u[7:end]
    # dX, dY, dC
    dXs = @view du[1:2:5]
    dYs = @view du[2:2:6]
    dC = @view du[7:end]

    #First moments
    dXs = -γ0/2*Xs + (ω - ωm)*Ys
    dYs = -γ0/2*Ys - (ω - ωm)*Xs

    # Covariance matrix
    for M = 1:2, N = 1:2
        for i in 1:3, j in 1:3
            Mi = M == 1 ?  2*(i-1)+1 : 2*i
            Nj = N == 1 ?  2*(j-1)+1 : 2*j
            Mconji = M == 1 ?  2*i : 2*(i-1)+1
            Nconjj = N == 1 ?  2*j : 2*(j-1)+1

            dC[getindex(Mi, Nj)] = - (γ0[i]+γ0[j])/2 * C[getindex(Mi,Nj)] + δ(Mi,Nj)*γ0*(nth+0.5) + δ(M,N)*sqrt(Γba[i]*Γba[j]) 
                     + (-1)^δ(M,2)*(ω[i] - ωm)*C[getindex(Mconji, Nj)] + (-1)^δ(N,2)*(ω[j] - ωm)*C[getindex(Mi, Nconjj)] 
                     -4*(sum([sqrt(Γmeas[k])*C[getindex(Mi, 2*k)] for k in 1:3]))*(sum([sqrt(Γmeas[l])*C[getindex(Nj, 2*k)] for Xl in 1:3]))
                     -4*(sum([sqrt(Γmeas[k])*C[getindex(Mi, 2*(k-1)+1)] for k in 1:3]))*(sum([sqrt(Γmeas[l])*C[getindex(Nj, 2*(k-1)+1)] for Yl in 1:3]))
        end
    end
    nothing
end

function moments_infogain(u,p,t)
    k, η = p[3:4]
    Vx, Cxp = @view u[3:2:5]
   sqrt(8*k*η)*SA[Vx, Cxp, 0, 0, 0]
end

function infogain_3modes(du,u,p,t)

end



# Parameters TODO: Convert to vectors
ωm = 2π*1.1
Q = 1e7
γ0 = ωm / Q

# Measurement "rate" TODO: Convert to vector
g0 = 2π*465
ncav = 1e5
g = g0*sqrt(ncav)
κ = 2π*45e6
Γqba = 4*g^2/κ /1e6
k = Γqba/2
η = 0.9

# Constants
ħ = 1.05457182e-34 * 1e9^2 / 1e6
kB = 1.380649e-23 * 1e9^2 / 1e6^2
T = 0#1e-3
const nth = 1/(exp(ħ*ωm/(kB*T))-1)
Vx_th = Vp_th = nth + 0.5

# Initial conditions and SDE parameters TODO: Change initial conditions to correct size
p = (γ0, ω, k, η)
u0 = SA[0, 0, Vx_th, Vp_th, 0.0]
tspan = (0.0, 100) # μs

prob= SDEProblem(moments_evolution, moments_infogain, u0, tspan, p)
@time sol = solve(prob,
    SOSRI(),
    saveat = 0.1/ω,
    dt = 1e-1,
    maxiters = tspan[2]*1e7,
    progress = true,
    save_idxs = 1
);

# Data visualization
pvars = plot(sol.t, sol[1,:], label = L"\langle x \rangle", xlabel = "t [μs]", ylabel = "Nondimensional")
plot!(sol.t, sol[2,:], label = L"\langle p \rangle")