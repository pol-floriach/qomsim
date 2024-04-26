# Exact numerical solution of the master equation of a damped harmonic oscillator with backaction

using LinearAlgebra, Kronecker, FastExpm, SparseArrays

# Functions for ladder operators
function anihilation(n)
    a = zeros(n,n)
    for i in 1:n-1
        a[i,i+1] = sqrt(i)
    end
    return a
end

function creation(n)
    at = zeros(n,n)
    for i in 1:n-1
        at[i+1,i] = sqrt(i)
    end
    return at
end

# ⊗(x,y) = kron(x,y)

# Parameters
ω_mech = 2π*1.1
Q = 1e7
γ₀ = ω_mech / Q

# Constants
ħ = 1.05457182e-34 * 1e9^2 / 1e6
kB = 1.380649e-23 * 1e9^2 / 1e6^2
T = 1e-3
n = 1/(exp(ħ*ω_mech/(kB*T))-1)
k = γ₀*n; k2 = 2*k

N = 100
# Operators
a = Bidiagonal(anihilation(N+1), :U)
at = Bidiagonal(creation(N+1),:L)
Nop = Diagonal(at*a)
H = Diagonal(ω_mech*Nop)
x = a+at
Id = I(N+1)

# Initial conditions

ρ0 = zeros(ComplexF64,N+1,N+1); ρ0[1,1] = 1;
ρ0vec = vec(ρ0)


# Time evolution: dₜρ̂ = 
O = -im*ω_mech*(Nop⊗Id - Id⊗Nop) + γ₀/2*Id⊗Id + γ₀*n*at⊗at + γ₀*(n+1)*a⊗a - γ₀/2*(1+2*n)*(Nop⊗Id + Id⊗Nop + Id⊗Id) + k2*(x⊗x - 0.5*(x^2⊗Id + Id⊗x^2))
Os = sparse(O)


ρtvec = similar(ρ0vec)
vecN = vec(at*a)# per fer Tr(A*B) = vec(A)†*vec(B)

# function timeevolution(tvec, Os, ρ0vec, )
it = 0
ns = Float64[]
for t in 1:1e6:1e7
    it += 1
    mul!(ρtvec, fastExpm(Os*t),ρ0vec)
    push!(ns,real(dot(ρtvec,vecN)))
    println("Iteration i = ", it)
end
# end





#heatmap(real.(reshape(ρtvec,(N+1,N+1))))

tvec = collect(1:1e6:1e7)
scatter!(tvec, ns, 
        xlabel = "t [μs]", 
        ylabel = "<n>", 
        label = "Exact numerical(N = 100), w/ba",
        ms = 2,
        msw = 0)

A = Nop ⊗ Id
B = Id⊗Nop

C = x⊗x - 0.5*(x^2⊗Id + Id⊗x^2)
