using ModelingToolkit, DifferentialEquations
@parameters t sig = 10 rho = 28.0 beta = 8 / 3
@variables x(t) = 100 y(t) = 1.0 z(t) = 1 u(t) = 1
D = Differential(t)

eqs = [D(x) ~ sig * (y - x),
    D(y) ~ x * (rho - z) - y,
    D(z) ~ x * y - beta * z,
    0 ~ u]

sys = ODESystem(eqs; tspan=(0, 100), name=:lorenz)

ssys = structural_simplify(sys)
prob = ODEProblem(sys)
prob2 = ODEProblem(ssys)
sol = solve(prob)
sol2 = solve(prob2)
fn = "lorenz_scimlstyle.jl"
write(fn, sys)
include(fn)

# function 