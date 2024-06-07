# Script for time evolution of the moments and covariances (SDE system) including 3 modes.TODO: Add appropriate initial conditions
using DifferentialEquations, Plots, LaTeXStrings

# Function to obtain O_k index of matrix (in a row or column)
function right_index(O,k,target= 1)
    return O == target ?  2k-1 : 2k
end
# Function to get index of a matrix stored in an array in column-major order (for 6x6 matrix)
colmaj(row, col) = row + (col-1)*6

# Kronecker delta function for better legibility
δ(x,y) = ==(x,y)

# Time evolution w/o information gain. TODO If too many allocations: try with staticarray / add Mi, Nj etc as parameters and mutate
function moments_evolution_3modes(du, u, p, t)
    # Parameters
    @views γ0, ωm, ωs, Γba ,Γmeas = p
    # X, Y, C
    Xs = @view u[1:2:5]
    Ys = @view u[2:2:6]
    C = @view u[7:end]
    # dX, dY, dC
    dXs = @view du[1:2:5]
    dYs = @view du[2:2:6]
    dC = @view du[7:end]

    #First moments
    @. dXs = -γ0/2*Xs + (ωs - ωm)*Ys
    @. dYs = -γ0/2*Ys - (ωs - ωm)*Xs

    # Covariance matrix
    for M = 1:2, N = 1:2
        for i in 1:3, j in 1:3
            # Getting the right indices
            Mi = right_index(M,i)
            Nj = right_index(N,j)
            Mconji = right_index(M,i,2)
            Nconjj = right_index(N,j,2)

            dC[colmaj(Mi, Nj)] = - (γ0[i]+γ0[j])/2 * C[colmaj(Mi,Nj)] + δ(Mi,Nj)*γ0[i]*(nth+0.5) + δ(M,N)*sqrt(Γba[i]*Γba[j]) +
                     + (-1)^δ(M,2)*(ωs[i] - ωm)*C[colmaj(Mconji, Nj)] + (-1)^δ(N,2)*(ωs[j] - ωm)*C[colmaj(Mi, Nconjj)] +
                     -4*(sum([sqrt(Γmeas[k])*C[colmaj(Mi, 2k-1)] for k in 1:3]))*(sum([sqrt(Γmeas[l])*C[colmaj(Nj, 2l-1)] for l in 1:3])) +
                     -4*(sum([sqrt(Γmeas[k])*C[colmaj(Mi, 2k)] for k in 1:3]))*(sum([sqrt(Γmeas[l])*C[colmaj(Nj, 2l)] for l in 1:3]))
        end
    end
    nothing
end

function infogain_3modes(du,u,p,t)
    C = @view u[7:end]
    
    du[1:2:5,1] = 2*[sum([sqrt(Γmeas[j])*C[colmaj(Xi,2*(j-1)+1)] for j in 1:3]) for Xi in 1:2:5] # dWx in d<Xi>
    du[1:2:5,2] = 2*[sum([sqrt(Γmeas[j])*C[colmaj(Xi,2*j)]       for j in 1:3]) for Xi in 1:2:5] # dWy in d<Xi>
    du[2:2:6,1] = 2*[sum([sqrt(Γmeas[j])*C[colmaj(Yi,2*(j-1)+1)] for j in 1:3]) for Yi in 2:2:6] # dWx in d<Yi>
    du[2:2:6,2] = 2*[sum([sqrt(Γmeas[j])*C[colmaj(Yi,2*j)]       for j in 1:3]) for Yi in 2:2:6] # dWy in d<Yi>
end



# Parameters
ωm = 2π*1.1
ωs = [2π*1.1, 2π*(1.1 - 1e-2), 2π*(1.1 + 1e-2)]
Qs  = [1e7, 1e7, 1e7]
γ0 = ωs ./ Qs


# Measurement "rates"
g0 = 2π*465
ncav = 1e5
g = g0*sqrt(ncav)
κ = 2π*45e6
Γqba = 4*g^2/κ /1e6
η = 0.9

Γba_vec = [Γqba, Γqba, Γqba]
Γmeas_vec = [Γqba*η, Γqba*η, Γqba*η]


# Constants
ħ = 1.05457182e-34 * 1e9^2 / 1e6
kB = 1.380649e-23 * 1e9^2 / 1e6^2
T = 0#1e-3
const nth = 1/(exp(ħ*ωm/(kB*T))-1)
Vx_th = Vp_th = nth + 0.5

# Initial conditions and SDE parameters TODO: Change initial conditions to correct size
p = (γ0, ωm, ωs, Γba_vec ,Γmeas_vec)

u0 = vcat(zeros(6), vec(sol_ss))

tspan = (0.0, 1000) # μs

prob= SDEProblem(moments_evolution, moments_infogain, u0, tspan, p)
@time sol = solve(prob,
    SOSRI(),
    saveat = 0.1/ω,
    dt = 1e-1,
    maxiters = tspan[2]*1e7,
    progress = true,
    save_idxs = 1
);

prob= ODEProblem(moments_evolution_3modes, u0, tspan, p)
@time sol = solve(prob,
    Tsit5(),
    saveat = 0.1/ωm,
    dt = 1e-1,
    maxiters = tspan[2]*1e7,
    progress = true,
);

# Data visualization
pvars = plot(sol.t, sol[1,:], label = L"\langle x \rangle", xlabel = "t [μs]", ylabel = "Nondimensional")
plot!(sol.t, sol[2,:], label = L"\langle p \rangle")


# Plotting time evolution of the variances (of the mode we aim to study, x1)

pcovs = plot(sol.t, sol[7,:], label = L"C_{X1X1}")
plot!(sol.t, sol[8,:], label = L"C_{X1Y1}")
`plot!(sol.t, sol[9,:], label = L"C_{X1X2}")
plot!(sol.t, sol[10,:], label = L"C_{X1Y2}")
plot!(sol.t, sol[11,:], label = L"C_{X1X3}")
plot!(sol.t, sol[12,:], label = L"C_{X1Y3}")