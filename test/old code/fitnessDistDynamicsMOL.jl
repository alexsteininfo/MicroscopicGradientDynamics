### title: numerical evalutation of PDEs with MethodsofLines
### author: Nathaniel Mon Pere

using ModelingToolkit, MethodOfLines, OrdinaryDiffEq, DomainSets

# Parameters, variables, and derivatives
@parameters t x

@variables n(..)
Dt = Differential(t)
Dx = Differential(x)
Dxx = Differential(x)^2

# parameters
w=0.1
b0 = 1.2
d0 = 0.2
K = 1000
ri = 0.
α = 10.
h = 1.

# birth and death rates
# birthRate(n,b0,K) = b0*(1 − n/K)
birthRate(n,b0,K) = b0
deathRate(h,u,d0,ri,α) = d0 + h/(ri + u*α)

# 1D PDE and boundary conditions
eq  = Dt(n(t, x)) ~ 
    (birthRate(n(t,x), b0, K) - deathRate(h, x, d0, ri, α))*n(t,x) + w/2 * Dxx(n(t, x))
bcs = [n(0, x) ~ 1.,
        Dx(n(t, 0)) ~ 0.0,
        Dx(n(t, 1)) ~ 0.0]

# Space and time domains
domains = [t ∈ (0.0, 10.0),
        x ∈ (0.0, 1.0)]

# PDE system
@named pdesys = PDESystem(eq, bcs, domains, [t, x], [n(t, x)])

# Method of lines discretization
dx = 0.01
# order = 2
discretization = MOLFiniteDifference([x => dx],t)

# Convert the PDE problem into an ODE problem
prob = discretize(pdesys, discretization)

# Solve ODE problem
using OrdinaryDiffEq
sol = solve(prob, Tsit5(), saveat=0.2)

##

using CairoMakie

# Plot results and compare with exact solution
discrete_x = sol[x]
discrete_t = sol[t]
soln = sol[n(t, x)]

fig1 = Figure()
Axis(
    fig1[1,1],
    xlabel="phenotype",
    ylabel="density of cells",
)
for t in [1,10,20,30,40,50]
    lines!(discrete_x, soln[t,:], label="time $t")
end
axislegend(position=:lt)
display(fig1)