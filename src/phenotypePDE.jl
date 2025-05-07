
using ModelingToolkit, MethodOfLines, OrdinaryDiffEq, DomainSets, Parameters

# abstract type Grid end

# struct UniformGrid <: Grid
#     u::Vector{Float64}
# end

# struct ChebyshevGrid <: Grid
#     L::Int
# end

# struct SymmetricExponentialGrid <: Grid
#     L::Int
#     ξ::Float64
# end

# abstract type GridType end

# struct UniformGridType <: Grid end

# struct ChebyshevGridType <: Grid end

# struct SymmetricExponentialGridType <: Grid end

function chebyGridpoint(i, L)
    return 1/2*(1 - cos(π*i/L))
end

function symExpoGridpoint(i, L, ξ)
    s = 2i/L-1
    1/2*( 1 + tanh(ξ*s)/tanh(ξ) )
end

function makeGrid(gridParams::NamedTuple)
    # gridParams = (type=..., L=..., varargs...)
    type, L = (gridParams.type, gridParams.L)
    _u = Vector{Float64}(undef, L)
    if type=="uniform"
        _u .= [i/(L-1) for i in 0:(L-1)]
    elseif type=="chebyshev"
        _u .= [chebyGridpoint(i, L-1) for i in 0:(L-1)]
    elseif type=="exponential"
        _u .= [symExpoGridpoint(i, L-1, gridParams.ξ) for i in 0:(L-1)]
    else
        error("Error: invalid grid type specified.")
    end
    return _u
end

function phenotype1DPDE(modelParams::NamedTuple, ctrlParams::NamedTuple)
    @unpack n0, b0, d0, l, ρUp, ρDown, K = modelParams
    @unpack t0, tF, gridParams, approxOrder, saveat = ctrlParams
    @parameters t u
    @variables n(..) totalPop(..)
    Dt = Differential(t)
    Du = Differential(u)
    Duu = Differential(u)^2
    uMin = 0.
    uMax = 1.
    Iu = Integral(u in DomainSets.ClosedInterval(uMin, uMax))

    uGrid = makeGrid(gridParams)
    
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
        # l^2/2*(ρ₊(uMin) + ρ₋(uMin))*Du(n(t,uMin)) - ( l*(ρ₊(uMin)-ρ₋(uMin)) - l^2/2*(Du(ρ₊(uMin))+Du(ρ₋(uMin))) ) * n(t,uMin) ~ 0.0,
        # l^2/2*(ρ₊(uMax) + ρ₋(uMax))*Du(n(t,uMax)) - ( l*(ρ₊(uMax)-ρ₋(uMax)) - l^2/2*(Du(ρ₊(uMax))+Du(ρ₋(uMax))) ) * n(t,uMax) ~ 0.0,
        l^2/2*(ρ₊(uMin) + ρ₋(uMin)) * Du(n(t,uMin)) - l*(ρ₊(uMin) - ρ₋(uMin))*n(t,uMin) ~ 0.0,
        l^2/2*(ρ₊(uMax) + ρ₋(uMax)) * Du(n(t,uMax)) - l*(ρ₊(uMax) - ρ₋(uMax))*n(t,uMax) ~ 0.0,
        Du(totalPop(t,uMin)) ~ 0,
        Du(totalPop(t,uMax)) ~ 0,
    ]

    domains = [
        t ∈ Interval(t0, tF),
        u ∈ Interval{:open,:open}(uMin, uMax),
    ]
    @named pdesys = PDESystem(equations, bcs, domains, [t, u], [n(t,u), totalPop(t,u)])
    discretization = MOLFiniteDifference([u => uGrid],t;approx_order=approxOrder,
    # grid_align=edge_align,
    should_transform=false,
    )
    prob = discretize(pdesys, discretization)
    solver = Tsit5()
    # solver = Rodas5P()
    # solver = Rodas4P()
    # solver = Trapezoid()
    sol = solve(prob, solver, saveat=saveat)
    _u = sol[u]
    _t = sol[t]
    n_t_u = sol[n(t, u)]
    return _t, _u, n_t_u
end