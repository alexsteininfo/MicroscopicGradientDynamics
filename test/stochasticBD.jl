# title: stochastic birth-death process with phenotypic switching
# author: alexander stein

using CairoMakie

######################
### Initialization ###
######################

### Initial condition
kspace = 1:100
n_init = [100 for k in kspace]
t_init = 0.0
t_eps = 1.0

### Stopping criteria
pop_max = 500000 #+Inf
pop_min = 0
t_max = 5.0

### Growth functions
b0 = 1.2
b(k) = b0

d0 = 0.2
l = 0.01
ri = 0.0
alpha = 10.0
h = 1.0
d(k) = d0 + h/(ri+k*l*alpha) #d0 + h/(ri + u*α)

### Transition functions
rp0 = 0.5
rp(k) = rp0

rm0 = 0.5
rm(k) = rm0


################
### Plotting ###
################

function plot_distribution()
    f = Figure()
    ax = Axis(f[1,1],
    #xticks = ([0,1000,2000,3000,4000]),
    xlabel = "Trait", xlabelsize = 30, xticklabelsize = 25,
    ylabel = "Density", ylabelsize = 30, yticklabelsize = 25
    )

    lines!(ax, kspace, nseries[1])
    lines!(ax, kspace, nseries[end])

    #hist!(exponents, color = (ColPaired[2], 0.6), normalization=:pdf, bins=20)
    #axislegend(ax, labelsize=25, framecolor = :white)

    save("MicroscopicGradientDynamics/test/figures/distributions.png", f)
end


#######################
### Run simulations ###
#######################

include("gillespiealgorithm.jl")

#run(maxtime::Float64 = Inf, minsize::Int64 = 0, maxsize::Int64 = Inf, maxsteps::Int64 = 1e6)
@time tseries, nseries = run(t_max, pop_min, pop_max, Int32(1e6))

plot_distribution()







