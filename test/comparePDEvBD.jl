include("../src/phenotypeSwitching.jl")
using .PhenotypeSwitching

## ----------------------------------------
#region - Solve PDE model

# these parameters control the model
modelParams = (
    b0 = 0.5, # birth rate
    d0 = 0.1, # death rate
    l = 0.02, #√(0.1), # discrete->continuous conversion parameter
    ρUp = 0.2,#0.1, # switch rate up
    ρDown = 0.5,#0.8, # switch rate down
    n0 = 2550., # initial population size
    K = 10000, # max population size
)
# these parameters control the solver
ctrlParams = (
    gridParams = (
        type="exponential",
        L=40,
        ξ=4,
    ),
    approxOrder=2,
    t0=0., # start time
    tF=10., # stop time
    saveat=1., # time step to save at
)

@time _t, _u, nPDE_t_u = PhenotypeSwitching.phenotype1DPDE(modelParams, ctrlParams)

#endregion


## ----------------------------------------
    #region - Solve Birth-Death model

# these parameters control the model
modelParams = (
    b0 = 0.5, # birth rate
    d0 = 0.1, # death rate
    ρUp = 0.2, # switch rate up
    ρDown = 0.5, # switch rate down
    u0 = [50 for i in 1:51], # initial population size, vector of length L
    K = 10000, # max population size
)
# these parameters control the solver
ctrlParams = (
    t0=0., # start time
    tF=10., # stop time
    saveat=1., # time step to save at
)

#@time _t1, _u1, n_t_u1, rn, N_t = phenotypeBD(modelParams, ctrlParams)
@time _t, _v, nBDMean_t_v = phenotypeBDMean(modelParams, ctrlParams, 100)


range(0,1,length=51)
#endregion

## ----------------------------------------
#region - Plotting

using CairoMakie

colorS = Makie.wong_colors()

fig1 = Figure()
Axis(
    fig1[1,1],
    xlabel="phenotype",
    ylabel="density of cells",
)
for t in [0,5,10]
    tInd = findfirst(_t.==t)
    lines!(_u, nPDE_t_u[tInd,:], label="PDE", color=colorS[1])
    lines!(_v, nBDMean_t_v[tInd,:], label="BD", color=colorS[2])
end
# ylims!(0, 2E4)
axislegend(position=:rt, merge=true)
# save("test/figures/PDEmodel.png",fig1)
display(fig1)

#endregion