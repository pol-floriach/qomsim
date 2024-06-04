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
function moments_evolution_3modes(u, p, t) # TODO: fix bug that says C or dC are float64s
    # Parameters
    γ0, ωm, ωs, Γba ,Γmeas = p
    # X, Y, C
    Xs = @view u[1:2:5]
    Ys = @view u[2:2:6]
    C = reshape(u[7:end],6,6)
    dC = similar(C)
    # dX, dY, dC
    #First moments
    dXs = @. -γ0/2*Xs + (ωs - ωm)*Ys
    dYs = @. -γ0/2*Ys - (ωs - ωm)*Xs

    # Covariance matrix
    for M = 1:2, N = 1:2
        for i in 1:3, j in 1:3
            # Getting the right indices
            Mi = right_index(M,i)
            Nj = right_index(N,j)
            Mconji = right_index(M,i,2)
            Nconjj = right_index(N,j,2)
            dC[Mi,Nj] = - (γ0[i]+γ0[j])/2 *  C[Mi,Nj] + δ(Mi,Nj)*γ0[i]*(nth+0.5) + δ(M,N)*sqrt(Γba[i]*Γba[j]) +
                     + (-1)^δ(M,2)*(ωs[i] - ωm)*C[Mconji, Nj] + (-1)^δ(N,2)*(ωs[j] - ωm)*C[Mi, Nconjj] +
                     - 4*(sum([sqrt(Γmeas[k])*C[Mi, 2k] for k in 1:3]))*(sum([sqrt(Γmeas[l])*C[Nj, 2l] for l in 1:3])) +
                     - 4*(sum([sqrt(Γmeas[k])*C[Mi, 2k-1+1] for k in 1:3]))*(sum([sqrt(Γmeas[l])*C[Nj, 2l-1] for l in 1:3]))
        end
    end
    return vcat(dXs, dYs, vec(dC))
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

# Initial conditions for finding the root
p = (γ0, ωm, ωs, Γba_vec ,Γmeas_vec)

u0_newton = zeros(42)

prob_ss= NonlinearProblem(moments_evolution_3modes, u0_newton, p)
@time sol = solve(prob_ss, NewtonRaphson(),
    reltol = 1e-10,
    abstol = 1e-10,
);
