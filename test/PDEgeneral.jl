# title: numerical evalutation of continuous trait space PDE
# author: Nathaniel Mon Pere

using ModelingToolkit, MethodOfLines, OrdinaryDiffEq, DomainSets, Parameters

# numerical parameters
# iginore the whole packing/unpacking thing for now. Becomes useful when everything gets more complex
params = (
    b0 = 0.5,
    d0 = 0.1,
    l = 0.01, #√(0.1),
    p = 0.0,
)
@unpack b0, d0, l, p = params
paramsAdaptiveTherapy = (
    K = 1000,
    ri = 0.,
    α = 10.,
    h = 1.
)
@unpack K, ri, α, h = paramsAdaptiveTherapy

# Parameters, variables, and derivatives
@parameters t u

@variables n(..)
Dt = Differential(t)
Du = Differential(u)
Duu = Differential(u)^2

β(u) = b0 + 0.1*u
δ(u) = d0 #d0 + h/(ri + u*α) # adaptive therapy model
ρ₊(u) = 0.5
ρ₋(u) = 0.5

# general phenotype switching equation
eq  = Dt(n(t, u)) ~ 
    ( β(u) - δ(u) - l*Du(ρ₊(u)) - l*Du(ρ₋(u)) + l^2/2*Duu(ρ₊(u)) + l^2/2*Duu(ρ₋(u)) ) * n(t,u) + 
    l*( -ρ₊(u) + ρ₋(u) + l*Du(ρ₊(u)) + l*Du(ρ₋(u)) ) * Du(n(t,u)) +
    l^2/2*( ρ₊(u) + ρ₋(u) ) * Duu(n(t,u))

bcs = [
    n(0, u) ~ 50.,
    Du(n(t, 0)) ~ 0.0,
    Du(n(t, 1)) ~ 0.0
]

# Space and time domains
domains = [
    t ∈ (0.0, 20.0),
    u ∈ (0.0, 1.0),
]

# PDE system
@named pdesys = PDESystem(eq, bcs, domains, [t, u], [n(t, u)])

# Method of lines discretization
du = 0.01
discretization = MOLFiniteDifference([u => du],t)
# Convert the PDE problem into an ODE problem
prob = discretize(pdesys, discretization)
# Solve ODE problem
using OrdinaryDiffEq
sol = solve(prob, Tsit5(), saveat=0.2)

##

using CairoMakie

# Plot results and compare with exact solution
_u = sol[u]
_t = sol[t]
n_t_u = sol[n(t, u)]

fig1 = Figure()
Axis(
    fig1[1,1],
    xlabel="phenotype",
    ylabel="density of cells",
)
for i in 1:5
    lines!(_u, n_t_u[findfirst(_t.==i*4.),:], label="t=$t")
end
# axislegend(position=:lt)
save("/Users/stein02/Desktop/github/MicroscopicGradientDynamics/test/figures/PDE_distributions.png", fig1)
#display(fig1)
