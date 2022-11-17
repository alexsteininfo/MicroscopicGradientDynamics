using ModelingToolkit, MethodOfLines, OrdinaryDiffEq, DomainSets

##

# Parameters, variables, and derivatives
@parameters t x
@parameters γ, w
@variables n(..)
Dt = Differential(t)
Dx = Differential(x)
Dxx = Differential(x)^2

# 1D PDE and boundary conditions
eq  = Dt(n(t, x)) ~ γ*x*n(t,x) + w/2 * Dxx(n(t, x))
bcs = [n(0, x) ~ 1.,
        Dx(n(t, 0)) ~ 0.0,
        Dx(n(t, 1)) ~ 0.0]

# Space and time domains
domains = [t ∈ Interval(0.0, 10.0),
        x ∈ Interval(0.0, 1.0)]

# PDE system
@named pdesys = PDESystem(eq, bcs, domains, [t, x],[n(t, x)], [γ=>0.1, w=>1.0])

dx = 0.01
discretization = MOLFiniteDifference([x => dx], t)

# Method of lines discretization
# Need a small dx here for accuracy
dx = 0.01
# order = 2
discretization = MOLFiniteDifference([x => dx],t)

# Convert the PDE problem into an ODE problem
prob = discretize(pdesys, discretization)

# Solve ODE problem
using OrdinaryDiffEq
sol = solve(prob, Tsit5(), saveat=0.2)


# Plot results and compare with exact solution
discrete_x = sol[x]
discrete_t = sol[t]

soln = sol[n(t, x)]

##

using CairoMakie

fig1 = Figure()
Axis(
    fig1[1,1],
    xlabel="phenotype",
    ylabel="density of cells",
)
lines!(discrete_x, soln[25,:])
display(fig1)