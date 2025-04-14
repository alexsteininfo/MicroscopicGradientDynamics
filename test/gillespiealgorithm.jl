# title: gillsepie algorithm
# author: alexander stein

#################
### libraries ###
#################

using Random
using StatsBase: sample, Weights

#####################
### Step 0: Input ###
#####################

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


### Derived parameters
N_init = sum(n_init)
num_types = length(kspace)
num_rates = 4*num_types

### Sanity check
if(pop_max <= N_init || pop_min >= N_init || t_max <= t_init)
    println("ERROR: Ill defined stopping criteria")
elseif( all(k -> b(k,N) < 0, kspace) || all(k -> d(k,N) < 0, kspace) )
    println("ERROR: Ill defined growth rate functions")
elseif( all(k -> rp(k) < 0, kspace) || all(k -> rm(k) < 0, kspace) )
    println("ERROR: Ill defined transition rate functions")
end

##############################
### Step 1: Random numbers ###
##############################

function gil_rand(gil_rates, gil_sum)
    # Random variable: Time tau
    tau = -log(rand())/gil_sum  # exponential with scale sum(a)
    # Random variable: Reaction rate mu
    mu = sample(Weights(gil_rates)) 

    return tau, mu
end

function create_rates!()    
    gilrates_b = []
    gilrates_d = []
    gilrates_rp = []
    gilrates_rm = []

    for k in kspace 
        push!(gilrates_b, b(k,N)*n[k])
        push!(gilrates_d, d(k,N)*n[k])
        push!(gilrates_rp, rp(k)*n[k])
        push!(gilrates_rm, rm(k)*n[k])
    end

    # We don't allow transitions out of bound
    gilrates_rp[end] = 0
    gilrates_rm[1] = 0

    global gil_rates = Real[gilrates_b; gilrates_d; gilrates_rp; gilrates_rm]
    global gil_sum = sum(gil_rates)
end

function update_rates!(k_update, eventtype)

    # update the rates of the k-th type in case of any event
    gil_rates[k_update]             = b(k_update,N)*n[k_update]
    gil_rates[k_update+num_types]   = d(k_update,N)*n[k_update]
    gil_rates[k_update+2*num_types] = rp(k_update)*n[k_update]
    gil_rates[k_update+3*num_types] = rm(k_update)*n[k_update]

    
    if(eventtype == 3 && k_update!=num_types)   # up transition
        gil_rates[k_update+1]               = b(k_update+1,N)*n[k_update+1]
        gil_rates[k_update+1+num_types]     = d(k_update+1,N)*n[k_update+1]
        gil_rates[k_update+1+2*num_types]   = rp(k_update+1)*n[k_update+1]
        gil_rates[k_update+1+3*num_types]   = rm(k_update+1)*n[k_update+1]
    elseif(eventtype == 4 && k_update!=1)   # down transition
        gil_rates[k_update-1]               = b(k_update-1,N)*n[k_update-1]
        gil_rates[k_update-1+num_types]     = d(k_update,N)*n[k_update-1]
        gil_rates[k_update-1+2*num_types]   = rp(k_update-1)*n[k_update-1]
        gil_rates[k_update-1+3*num_types]   = rm(k_update-1)*n[k_update-1]
    end

    gil_rates[3*num_types] = 0
    gil_rates[3*num_types+1] = 0

    gil_sum = sum(gil_rates)

end

#####################################
### Step 2: Update time and state ###
#####################################

function update_type(mu)
    if(mu <= num_types)
        k_update = mu
        eventtype = 1
    elseif(mu <= 2*num_types)
        k_update = mu-num_types
        eventtype = 2
    elseif(mu <= 3*num_types)
        k_update = mu-2*num_types
        eventtype = 3
    elseif(mu <= 4*num_types)
        k_update = mu-3*num_types
        eventtype = 4
    else
        println("ERROR: Ill defined gillespie update rate!")
    end

    return k_update, eventtype
end

function update_state!(tau, k_update, eventtype)
    # Update time t
    global t += tau

    # Update population sizes n
    if(eventtype == 1)       # birth
        n[k_update] += 1
        global N += 1
    elseif(eventtype == 2)   # death
        n[k_update] -= 1
        global N -= 1
    elseif(eventtype == 3)   # up transition
        n[k_update] -= 1
        n[k_update+1] += 1
    elseif(eventtype == 4)   # down transition
        n[k_update] -= 1
        n[k_update-1] += 1
    else
        println("ERROR: Ill defined update rules!")
    end
end



#############################
### Step 3: Read and save ###
#############################

function save!()
    push!(nseries, copy(n))
    push!(tseries, t)
end

#################################
### Step 4: Stopping criteria ###
#################################

function stop(maxtime, minsize, maxsize, maxsteps)
    stopping = false 
    if(N==0)
        println("POPULATION EXTINCT")
        stopping = true
    elseif(t>=maxtime)
        println("REACHED MAXIMUM TIME")
        stopping = true
    elseif(N<=pop_min)
        println("MINIMUM SIZE REACHED")
        stopping = true
    elseif(N>=pop_max)
        println("MAXIMUM SIZE REACHED")
        stopping = true
    elseif(i==maxsteps)
        println("REACHED MAXIMUM STEPS")
        stopping = true
    end

    return stopping
end


######################
### Run simulation ###
######################

function run(maxtime, minsize, maxsize, maxsteps)
    # Initialize state variables
    global n = copy(n_init)
    global N = N_init
    global t = t_init

    # Initialize saving
    global nseries = [copy(n)]
    global tseries = [t]

    # Initialize helping variables
    global t_save = t_init
    create_rates!()

    # Run through some steps
    global i = 1
    while i<=maxsteps
        global i += 1
        # Generate random numbers
        tau, mu = gil_rand(gil_rates, gil_sum)
        # Determine event type
        k_update, eventtype = update_type(mu)
        # Update time and population sizes
        update_state!(tau, k_update, eventtype)
        # Update gillsepie rates
        #update_rates!(gil_rates, gil_sum)
        create_rates!()
        #update_rates!(k_update, eventtype)
        # Test stopping criteria
        if(stop(maxtime, minsize, maxsize, maxsteps))
            save!()
            break
        end
        # Save after time t_eps
        if(t > t_save + t_eps)
            save!()
            t_save += t_eps
        end
    end
    return tseries, nseries
end


