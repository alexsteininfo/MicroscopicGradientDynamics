# title: Birth-death process with transitions
# author: alexander stein

### Initial condition
# kspace - allowed values of trait k
# n_init - initial pop sizes in state k 
# t_init - initial time for simulation
# t_eps - time steps at which pop size is saved

### Stopping criteria
# pop_max
# pop_min
# t_max

### Growth functions
# b(k,N) birth function
# d(k,N) death function

### Transition functions
# rp(k) upward transition function
#   implementation forces rp(k_max)=0
# rm(k) downward transition function
#   implementation forces rm(k_min)=0

using DifferentialEquations, Catalyst

function phenotypeBD(modelParams::NamedTuple, ctrlParams::NamedTuple)
    
    @unpack u0, b0, d0, ρUp, ρDown, K = modelParams
    @unpack t0, tF, saveat = ctrlParams

    L = length(u0)
    _u = range(0.0, 1.0, L)

    t = default_t()
    @species (S(t)[1:L]) N(t) E(t)
    @parameters birthrates[1:L] deathrates[1:L] competitionrates[1:L] rhominus[1:L] rhoplus[1:L] zerorate

    # Empty reactions to fix the order of reactants in the solution
    reactions = Reaction[]
    for i in 1:L
        push!(reactions, Reaction(zerorate, [S[i]], [S[i]], [1], [2]))
    end
    push!(reactions, Reaction(zerorate, [N], [N], [1], [2]))
    push!(reactions, Reaction(zerorate, [E], [E], [1], [2]))
    # Actual reactions
    for i in 1:L
        push!(reactions, Reaction(birthrates[i], [S[i]], [S[i],N], [1], [2,1]))
        push!(reactions, Reaction(deathrates[i], [S[i]], [S[i],N], [1], [0,-1]))
        #push!(reactions, Reaction(competitionrates[i], [S[i],N], [S[i],N], [1], [0,-1]))
        push!(reactions, Reaction(competitionrates[i], [S[i],N], [S[i],N], [1,1], [0,0]))
        #push!(reactions, Reaction(treatmentrates[i], [S[i]], [S[i],N], [1], [0,-1]))
        #i!=1 ? push!(reactions, Reaction(rhominus[i], [S[i]], [S[i-1]])) : push!(reactions, Reaction(rhoplus[L], [E], [E], [1], [2]))
        #i!=L ? push!(reactions, Reaction(rhoplus[i], [S[i]], [S[i+1]])) : push!(reactions, Reaction(rhominus[1], [E], [E], [1], [2]))
    end

    for i in 2:L
        push!(reactions, Reaction(rhominus[i], [S[i]], [S[i-1]]))
    end

    for i in 1:(L-1)
        push!(reactions, Reaction(rhoplus[i], [S[i]], [S[i+1]]))
    end

    @named rn = ReactionSystem(reactions, t)
    rn = complete(rn)

    p = (:birthrates => [b0+0.1*_u[i] for i in 1:L],
        :deathrates => [d0 for i in 1:L],
        :rhominus => [ρDown for i in 1:L],
        :rhoplus => [ρUp for i in 1:L],
        :competitionrates => [(b0+0.1*_u[i]-d0)/K for i in 1:L],
        #:competitionrates => [0.0 for i in 1:L],
        :zerorate => 0.0)

    u0 = [u0; sum(u0); 1]
    tspan = (t0, tF)

    prob = DiscreteProblem(rn, u0, tspan, p)

    jump_prob = JumpProblem(rn, prob, Direct(), save_positions = (false,false))

    sol = solve(jump_prob, SSAStepper(); saveat = saveat)
    _t = sol[t]
    n_t_u = reduce(vcat, (v' for v in sol[S])).*L

    N_t = sol[N]

    return _t, _u, n_t_u, rn, N_t
end


function phenotypeBD_mean(modelParams::NamedTuple, ctrlParams::NamedTuple, realisations)
    matrices = []
    
    _t, _u, n_t_u, rn, N_t = phenotypeBD(modelParams, ctrlParams)
    push!(matrices, n_t_u)
    for i in 2:realisations
        _t, _u, n_t_u, rn, N_t = phenotypeBD(modelParams, ctrlParams)
        push!(matrices, n_t_u)
    end

    sum_matrix = matrices[1]
    for i in 2:realisations
        sum_matrix += matrices[i]
    end

    n_t_u_mean = sum_matrix/realisations

    return _t, _u, n_t_u_mean
end

## ----------------------------------------
#region - Solve model

# these parameters control the model
modelParams = (
    b0 = 0.5, # birth rate
    d0 = 0.1, # death rate
    #l = 0.01, #√(0.1), # discrete->continuous conversion parameter
    ρUp = 0.01,#0.1, # switch rate up
    ρDown = 0.08,#0.8, # switch rate down
    #n0 = 50., # initial population size
    u0 = [50 for i in 1:51], # initial population size, vector of length L
    K = 10000, # max population size
)
# these parameters control the solver
ctrlParams = (
    #du=0.05, # discretization in the phenotype space
    t0=0., # start time
    tF=10., # stop time
    saveat=1., # time step to save at
)

#@time _t1, _u1, n_t_u1, rn, N_t = phenotypeBD(modelParams, ctrlParams)

@time _t, _u, n_t_u_mean = phenotypeBD_mean(modelParams, ctrlParams, 100)

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
    #lines!(_u1, n_t_u1[tInd,:], label="t=$t")
    lines!(_u, n_t_u_mean[tInd,:], label="t=$t")
end
axislegend(position=:lt)
#display(fig1)
save("test/figures/BDmodel.png",fig1)


#endregion