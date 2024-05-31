using DifferentialEquations, StaticArrays, Plots, LaTeXStrings

function moments_evolution_multimode(u,p,t)
    @views γ0, ω, k, η = p
    @views Xi, Yi, VXi, VYi, CXiYi = u
    #First moments
    dXi    = -γ0/2*Xi + (ωi - ωm)*Yi
    dYi    = -γ0/2*Yi - (ωi - ωm)*Xi
    # Covariance matrix

    # Precompute sum terms TODO: Check if we are counting XiXi, YiYi, XiYi. Also change CYiXj to CXjYi (bc its symmetrized covariance).
    sum_CXiXj = sum([sqrt(Γ_ba[j]*η_j)*CXiXj[j] for j in eachindex(CXiXj)])
    sum_CXiYj = sum([sqrt(Γ_ba[j]*η_j)*CXiYj[j] for j in eachindex(CXiYj)])
    sum_CYiXj = sum([sqrt(Γ_ba[j]*η_j)*CYiXj[j] for j in eachindex(CYiXj)])
    sum_CYiYj = sum([sqrt(Γ_ba[j]*η_j)*CYiYj[j] for j in eachindex(CYiYj)])

    dVXi   = -γ0*VXi + 2*(ωi - ωm)*CXiYi + γ0*(nth+0.5) + Γba_i - 4*((sqrt(Γba_i*η_i)*VXi + sum_CXiXj)^2 + (sqrt(Γba_i)*CXiYi + sum_CXiYj)^2)
    dVYi   = -γ0*VYi - 2*(ωi - ωm)*CXiYi + γ0*(nth+0.5) + Γba_i - 4*((sqrt(Γba_i*η_i)*VYi + sum_CYiYj)^2 + (sqrt(Γba_i)*CXiYi + sum_CYiXj)^2)
    dCXiYi  = ω*Vp - ω*Vx - γ0*Cxp - 8*k*η*Vx*Cxp
   
    # Crossed terms (seria un vector aixo)
   dCXiXj = 5
   dCXiYj = 6

    SA[dXi, dYi, dVXi, dVYi, dCXiYi, dCXiXj, dCXiYj]
end

function moments_evolution_3modes(u,p,t)
    @views γ0, ω, k, η = p
    @views X1, Y1, VX1, VY2, CX1Y1, CX1X2, CX1X3, CX1Y2, CX1Y3 = u
    #First moments
    dX1    = -γ0/2*X1 + (ω1 - ω1)*Y1
    dY1    = -γ0/2*Y1 - (ω1 - ω1)*X1
    
    # Covariance matrix
    dVX1   = -γ0*VX1 + 2*(ω1 - ω1)*CX1Y1 + γ0*(nth+0.5) + Γba_1 - 4*((sqrt(Γba_1*η_1)*VX1 + sqrt(Γba_2*η_2)*CX1X2 + sqrt(Γba_3*η_3)*CX1X3)^2 + (sqrt(Γba_i)*CX1Y1 + sqrt(Γba_2*η_2)*CX1Y2 + sqrt(Γba_3*η_3)*CX1Y3)^2)
    dVY1   = -γ0*VY1 - 2*(ω1 - ω1)*CX1Y1 + γ0*(nth+0.5) + Γba_1 - 4*((sqrt(Γba_1*η_1)*VY1 + sqrt(Γba_2*η_2)*CY1Y2 + sqrt(Γba_3*η_3)*CY1Y3)^2 + (sqrt(Γba_i)*CX1Y1 + sqrt(Γba_2*η_2)*CY1X2 + sqrt(Γba_3*η_3)*CY1X3)^2)
    dCX1Y1  = (ω1 - ω1)*(VY1 - VX1) - γ0*CX1Y1 - 4*((sqrt(Γba_1*η_1)*VX1 + sqrt(Γba_2*η_2)*CX1X2 + sqrt(Γba_3*η_3)*CX1X3)*(sqrt(Γba_1*η_1)*CY1X1 + sqrt(Γba_2*η_2)*CY1X2 + sqrt(Γba_3*η_3)*CY1X3))
                                               - 4*((sqrt(Γba_1*η_1)*VY1 + sqrt(Γba_2*η_2)*CY1Y2 + sqrt(Γba_3*η_3)*CY1Y3)*(sqrt(Γba_1*η_1)*CX1Y1 + sqrt(Γba_2*η_2)*CX1Y2 + sqrt(Γba_3*η_3)*CX1Y3))
   
    # Crossed terms (seria un vector aixo)
    dCX1X2 = 6
    dCX1X3 = 7
    dCX1Y2 = 8
    dCX1Y3 = 9
    dCY1Y2 = 10
    dCY1Y3 = 11

    SA[dXi, dYi, dVXi, dVYi, dCXiYi, dCXiXj, dCXiYj]
end



# Parameters
ωm = 2π*1.1
Q = 1e7
γ0 = ωm / Q

# Measurement "rate"
g0 = 2π*465
ncav = 1e5
g = g0*sqrt(ncav)
κ = 2π*45e6
Γqba = 4*g^2/κ /1e6
k = Γqba/2
η = 0.9
G = 1

# Constants
ħ = 1.05457182e-34 * 1e9^2 / 1e6
kB = 1.380649e-23 * 1e9^2 / 1e6^2
T = 0#1e-3
const nth = 1/(exp(ħ*ωm/(kB*T))-1)
Vx_th = Vp_th = nth + 0.5

# Initial conditions and SDE parameters
p = (γ0, ω, k, η)
u0 = SA[0, 0, Vx_th, Vp_th, 0.0]
tspan = (0.0, 100) # μs