using ModelingToolkit, DifferentialEquations
using SciMLBase, DiffEqBase, OrdinaryDiffEq, DifferentialEquations, ModelingToolkit, Symbolics
using Plots

@parameters t
@parameters p=-1.
@variables x(t) = 1 y(t) = 2

@variables x(t) = 1 y(..) = 2
D = Differential(t)

# eqs = [D(x) ~ p*x]
# @named sys = ODESystem(eqs)

# prob = ODEProblem(sys, [], (0, 1000))
# sol = solve(prob)
# plot(sol)

# eqs = [D(x) ~ p*x, 0 ~ 0]
eqs = [0 ~ 0]
eqs = [D(x) ~ p*x, y ~ 0]
@named sys = ODESystem(eqs, t)
ssys = structural_simplify(sys)
prob = ODEProblem(ssys, [], (0, 1000))
sol = solve(prob)
plot(sol)
