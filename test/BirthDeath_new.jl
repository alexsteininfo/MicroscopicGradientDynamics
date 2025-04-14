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

n = 5
t = default_t()
@species (S(t)[1:n]) N(t) E(t)
@parameters birthrates[1:n] deathrates[1:n] rhominus[1:n] rhoplus[1:n] zerorate

reactions = Reaction[]


for i in 1:n
    push!(reactions, Reaction(zerorate, [S[i]], [S[i]], [1], [2]))
end
push!(reactions, Reaction(zerorate, [N], [N], [1], [2]))
push!(reactions, Reaction(zerorate, [E], [E], [1], [2]))

for i in 1:n
    push!(reactions, Reaction(birthrates[i], [S[i]], [S[i],N], [1], [2,1]))
    push!(reactions, Reaction(deathrates[i], [S[i]], [S[i],N], [1], [0,-1]))
    push!(reactions, Reaction(competitionrates[i], [S[i],N], [S[i],N], [1], [0,-1]))
    #push!(reactions, Reaction(treatmentrates[i], [S[i]], [S[i],N], [1], [0,-1]))
    i!=1 ? push!(reactions, Reaction(rhominus[i], [S[i]], [S[i-1]])) : push!(reactions, Reaction(rhoplus[n], [E], [E], [1], [2]))
    i!=n ? push!(reactions, Reaction(rhoplus[i], [S[i]], [S[i+1]])) : push!(reactions, Reaction(rhominus[1], [E], [E], [1], [2]))
end

@named rn = ReactionSystem(reactions, t)
rn = complete(rn)

p = (:birthrates => [1.0 for i in 1:n], :deathrates => [0.1 for i in 1:n], :rhominus => [0.0 for i in 1:n], :rhoplus => [0.0 for i in 1:n], :zerorate => 0.0)
#p = (:birthrates => [1.0 for i in 1:n], :deathrates => [0.0 for i in 1:n])
u0 = [50 for k in 1:n]
u0 = [u0; sum(u0); 1]
tspan = (0, 5.0)

prob = DiscreteProblem(rn, u0, tspan, p)

jump_prob = JumpProblem(rn, prob, Direct(), save_positions = (false,false))

sol = solve(jump_prob, SSAStepper(); saveat = 1.0)
# specific timepoints can be accessed via sol(t), e.g. sol(1.5)

