# title: transient dynamics of the G-function approach
# author: alexander stein

### Import libraries
using DifferentialEquations # Solving differential equations
using Plots

### Parameter values
b0 = 0.45
d0 = 0.01
K = 1e4

m = 0.5
k = 2
b = 10

sigma = 0.1

param = [b0, d0, K, m, k, b, sigma]

# Initial condition
N0 = 100
u0 = 0

# Time interval
timeinterval = [0.0, 200.0]

initial_values = [N0, u0]

# G-function and its gradient
G(N,u) = b0*(1-N/K)-d0-m/(k+b*u)
dGdu(N,u) = b*m/(k+b*u)^2

### Define the ode system
function gfunction_ode!(ds, s, p, t)

    # unpack variables and parameters
    N, u = s
    b0, d0, K, m, k, b, sigma = p

    # define differential equation
    dN = N*G(N,u)
    du = sigma*dGdu(N,u)

    ds .= (dN, du)
end

### Solve ODEProblem
problem = ODEProblem(gfunction_ode!, initial_values, timeinterval, param)
solution = solve(problem, saveat = 0.1) 


function plot_popsize(solution)
    plot(solution.t, solution[1,:], label="x for lotka-volterra", color="red", linestyle=:solid, linewidth = 2)
    #plot!(yaxis=(:log, [y0/x0, :auto]), legend=:bottomright)
    plot!(xlabel="Time")
    plot!(ylabel="Population size (x)")
    plot!(legend=nothing)
    savefig("popsize_gfunction")
end

function plot_trait(solution)
    plot(solution.t, solution[2,:], lxxbel="y for lotka-volterra", color="blue", linestyle=:solid, linewidth = 2)
    #plot!(yaxis=(:log, [y0/x0, :auto]), legend=:bottomright)
    plot!(xlabel="Time")
    plot!(ylabel="Mean trait value (u)")
    plot!(legend=nothing)
    savefig("trait_gfunction")
end

plot_popsize(solution)
plot_trait(solution)





