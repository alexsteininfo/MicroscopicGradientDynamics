# title: birth-death process with phenotypic switching
# author: alexander stein

#################
### libraries ###
#################

using Random
using StatsBase: sample, Weights         # for sampling with Weights
using Plots
using Roots             # for fitting distribution
using DataStructures: counter


# Set random number seed
Random.seed!(1)


##############################
### Step 0: Initialization ###
##############################

# Initial population size
num_types = 100
X0 = [100 for i in 1:num_types]    # 10 cells of each of the 500 types at the beginning
u0 = [1.0*i/num_types for i in 1:num_types]

# Parameters
b0 = 0.45
d0 = 0.01
K = 1e4

m = 0.5
k = 2
b = 10

p = 0.6     # probability to change the state

birth(N) = b0*(1-N/K)
death(m,u) = d0 + m/(k+b*u)



##############################
### Step 1: Random numbers ###
##############################

function gil_rand(X)
    # Create a vector "a" with all possible reaction rates
    a = Real[]
    N = sum(X)
    
    for Xi in X     # Add birth rates
        a_birth = birth(N).*Xi
        push!(a, a_birth)
    end
    for (Xi, ui) in zip(X, u0)     # Add death rates
        a_death = death(m,ui).*Xi
        push!(a, a_death)
    end
    #println(a)

    # Random variable: Time tau
    tau = -log(rand())/sum(a)  # exponential with scale sum(a)

    # Random variable: Reaction rate mu
    indices = 1:size(a)[1]
    mu = sample(indices, Weights(a))

    return tau, mu
end


#########################################
### Step 2: Update time and reactants ###
#########################################

function update(X, t, tau, mu)
    # Update time t
    t += tau

    # Update population sizes
    if mu==1                 # it's a birth event on type 1
        q = rand()
        if q <= p/2
            X[mu+1] += 1
        else
            X[mu] += 1
        end
    elseif mu<=num_types-1  # it's a birth event on type mu
        q = rand()
        if q <= p/2
            X[mu+1]+=1
        elseif q <= p
            X[mu-1]+=1
        else
            X[mu] +=1
        end
    elseif mu==num_types    # it's a birth event on type num_types
        q = rand()
        if q <= p/2
            X[mu-1]+=1
        else
            X[mu] +=1
        end
    elseif mu<=2*num_types  # it's a death event of type mu
        X[mu-num_types] -= 1
    else
        println("Something wrong with the update rules!")
    end
    
    #println(mu)

    return X, t
end



#############################
### Step 3: Read and save ###
#############################

function save(Xseries, tseries, X, t)
    # Save X in Xseries and t in tseries after certain time steps delta
    if(t-tseries[1][end]>0.01) # Save evere 0.01 time units
        for i in 1:num_types
            push!(Xseries[i],X[i])
        end
        push!(tseries,t)
    end

    return Xseries, tseries
end


######################
### Run simulation ###
######################

function run(maxtime, maxsize, maxsteps)
    # Set initial values
    X = copy(X0)
    t = 0.0

    # Save population sizes and time
    Xseries = [[Xi] for Xi in X0]
    tseries = [t]

    # Run through some steps
    i = 1
    while i<=maxsteps
        # Generate random numbers
        tau, mu = gil_rand(X)
        # Update time and population sizes
        X, t = update(X, t, tau, mu)
        # Save the interesting data
        Xseries, tseries = save(Xseries, tseries, X, t)
        
        # Stopping when extinct
        if(sum(X)==0)
            println("POPULATION EXTINCT")
            break
        end
        # Stopping when reaching maxtime
        if(t>=maxtime)
            println("REACHED MAXIMUM TIME")
            break
        end
        # Stopping when reaching size zero
        if(i==maxsteps)
            println("REACHED MAXIMUM STEPS")
            break
        end

    end
    return tseries, Xseries
end


################
### Analysis ###
################

function create_Xend(Xseries)
    Xend = []
    for Xi in Xseries
        push!(Xend, Xi[end])
    end
    return Xend
end


################
### Plotting ###
################

function plot_distribution(dist)
    xaxis = u0
    yaxis = dist
    plt = plot(xaxis, yaxis)
    plot!(xlabel="Level of resistance u", ylabel="Abundance")
    savefig(plt, "dist")
    return plt
end

function plot_popsizes(tseries, Xseries)
    
    plt = plot()
    for Xi in Xseries
        plot!(tseries, Xi)
    end
    plot!(xlabel="Time", ylabel="Population size")
    plot!(legend=nothing)
    savefig(plt, "popsizes")
    return plt

end

############
### Main ###
############

maxtime = 5
maxsize = 10e7
maxsteps = 1e4
tseries, Xseries = run(maxtime, maxsize, maxsteps)


cXend = create_Xend(Xseries)
plt_dist = plot_distribution(cXend)
plt_sizes = plot_popsizes(tseries, Xseries)


