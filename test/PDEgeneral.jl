# title: numerical evalutation of continuous trait space PDE
# author: Nathaniel Mon Pere

using ModelingToolkit, MethodOfLines, OrdinaryDiffEq, DomainSets, Parameters

function phenotype1DPDE(modelParams::NamedTuple, ctrlParams::NamedTuple)
    @unpack n0, b0, d0, l, ρUp, ρDown, K = modelParams
    @unpack t0, tF, du, saveat = ctrlParams
    @parameters t u
    @variables n(..) totalPop(..)
    Dt = Differential(t)
    Du = Differential(u)
    Duu = Differential(u)^2
    uMin = 0.
    uMax = 1.
    Iu = Integral(u in DomainSets.ClosedInterval(uMin, uMax))
    
    β(u) = b0 + 0.1*u
    δ(u) = d0
    # δ(u) = d0 + h/(ri + u*α) # adaptive therapy model: not implemented here
    γ(u) = (β(u) - δ(u))*(1-totalPop(t,u)/K)
    ρ₊(u) = ρUp
    ρ₋(u) = ρDown
    equations  = [
        totalPop(t,u) ~ Iu(n(t,u)),
        Dt(n(t, u)) ~ 
        ( γ(u) - l*Du(ρ₊(u)) + l*Du(ρ₋(u)) + l^2/2*Duu(ρ₊(u)) + l^2/2*Duu(ρ₋(u)) ) * n(t,u) + 
        l*( -ρ₊(u) + ρ₋(u) + l*Du(ρ₊(u)) + l*Du(ρ₋(u)) ) * Du(n(t,u)) +
        l^2/2*( ρ₊(u) + ρ₋(u) ) * Duu(n(t,u))
    ]
    bcs = [
        n(0.0, u) ~ n0,
        Du(n(t, uMin)) ~ 0.0,
        Du(n(t, uMax)) ~ 0.0
    ]
    domains = [
        t ∈ Interval(t0, tF),
        u ∈ Interval(uMin, uMax),
    ]
    @named pdesys = PDESystem(equations, bcs, domains, [t, u], [n(t,u), totalPop(t,u)])
    # @time @named pdesys = PDESystem(equations, bcs, domains, [t, u], [n(t,u), totalPop(t,u)], [ρP, ρM], defaults=Dict([ρUp=>0.1, ρDown=>0.8]))
    discretization = MOLFiniteDifference([u => du],t)
    prob = discretize(pdesys, discretization)
    sol = solve(prob, Tsit5(), saveat=saveat)
    _u = sol[u]
    _t = sol[t]
    n_t_u = sol[n(t, u)]
    return _t, _u, n_t_u
end

## ----------------------------------------
#region - Solve model

# these parameters control the model
modelParams = (
    b0 = 0.5, # birth rate
    d0 = 0.1, # death rate
    l = 0.02, #√(0.1), # discrete->continuous conversion parameter
    ρUp = 0.1,#0.1, # switch rate up
    ρDown = 0.8,#0.8, # switch rate down
    n0 = 2550., # initial population size
    K = 10000, # max population size
)
# these parameters control the solver
ctrlParams = (
    du=0.05, # discretization in the phenotype space
    t0=0., # start time
    tF=10., # stop time
    saveat=1., # time step to save at
)

@time _t, _u, n_t_u = phenotype1DPDE(modelParams, ctrlParams)

#endregion

## ----------------------------------------
#region - Plot results
using CairoMakie

fig1 = Figure()
Axis(
    fig1[1,1],
    xlabel="phenotype",
    ylabel="density of cells",
)
for t in [0,5,10]
    tInd = findfirst(_t.==t)
    lines!(_u, n_t_u[tInd,:], label="t=$t")
end
axislegend(position=:lt)
save("test/figures/PDEmodel.png",fig1)
#display(fig1)

#endregion