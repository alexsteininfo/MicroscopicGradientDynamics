# title: stochastic birth-death process with phenotypic switching
# author: alexander stein

using CairoMakie

######################
### Initialization ###
######################

### Initial condition
kspace = 1:101
n_init = [50 for k in kspace]
#n_init = [0 for k in kspace]
#n_init[1] = 100
t_init = 0.0
t_eps = 5.0

### Stopping criteria
pop_max = 1000000 #+Inf
pop_min = 0
t_max = 20.0

### Growth functions
b0 = 0.5
b(k,N) = max(b0 + 0.1*l*k - N/10000, 0) #+ 0.01*k/100 #max(b0 - N/1000000.0,0)

d0 = 0.1
l = 1/length(kspace)
ri = 2.0
alpha = 10.0
h = 1.0
d(k,N) = d0 #+ h/(ri+k*l*alpha) #d0 + h/(ri + u*α)

### Transition functions
rp0 = 0.5
rp(k) = rp0

rm0 = 0.5
rm(k) = rm0


################
### Plotting ###
################

function plot_distribution(name, timepoint)
    f = Figure()
    ax = Axis(f[1,1],
    #xticks = ([0,1000,2000,3000,4000]),
    xlabel = "Phenotype, 𝑘", xlabelsize = 30, xticklabelsize = 25,
    ylabel = "Density, 𝑁ₖ", ylabelsize = 30, yticklabelsize = 25
    )

    #barplot!(ax, kspace, nseries[1])
    #barplot!(ax, kspace, nseries[2])
    #barplot!(ax, kspace, nseries[3])
    barplot!(ax, kspace, nseries[timepoint])
    ylims!(0,100)

    #hist!(exponents, color = (ColPaired[2], 0.6), normalization=:pdf, bins=20)
    #axislegend(ax, labelsize=25, framecolor = :white)

    save(name, f)
end


#######################
### Run simulations ###
#######################

include("gillespiealgorithm.jl")

#run(maxtime::Float64 = Inf, minsize::Int64 = 0, maxsize::Int64 = Inf, maxsteps::Int64 = 1e6)
@time tseries, nseries = run(t_max, pop_min, pop_max, Int32(1e6))
name = "test/figures/distribution1.png"
plot_distribution(name, 1)
name = "test/figures/distribution2.png"
plot_distribution(name, 2)
name = "test/figures/distribution3.png"
plot_distribution(name, 3)
#name = "test/figures/distribution4.png"
#plot_distribution(name, 4)
#name = "test/figures/distribution5.png"
#plot_distribution(name, 5)
#name = "test/figures/distribution6.png"
#plot_distribution(name, 6)
#name = "test/figures/distribution7.png"
#plot_distribution(name, 7)


#=
finalsize1 = Int32[]
for i in 1:50
    tseries, nseries = run(t_max, pop_min, pop_max, Int32(1e6))
    push!(finalsize1, sum(nseries[end]))
end

#println(sum.(nseries))


rp0 = 0.5
rp(k) = rp0
rm0 = 0.5
rm(k) = rm0

@time tseries, nseries = run(t_max, pop_min, pop_max, Int32(1e6))
name = "test/figures/distribution2.png"
plot_distribution()

finalsize2 = Int32[]
for i in 1:50
    tseries, nseries = run(t_max, pop_min, pop_max, Int32(1e6))
    push!(finalsize2, sum(nseries[end]))
end

println(sum.(nseries))

rp0 = 1.0
rp(k) = rp0
rm0 = 1.0
rm(k) = rm0

@time tseries, nseries = run(t_max, pop_min, pop_max, Int32(1e6))
name = "test/figures/distribution3.png"
plot_distribution()

println(sum.(nseries))

finalsize3 = Int32[]
for i in 1:50
    tseries, nseries = run(t_max, pop_min, pop_max, Int32(1e6))
    push!(finalsize3, sum(nseries[end]))
end

#hist(f[1, 1], data, bins = 10)
=#


