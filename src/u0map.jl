using ModelingToolkit, OrdinaryDiffEq
@parameters σ
@parameters t sig = 10 rho = 28.0 beta = 8 / 3
@variables x(t) = 100 y(t) = 1.0 z(t) = 1
D = Differential(t)

eqs = [D(x) ~ sig * (y - x),
    D(y) ~ x * (rho - z) - y,
    D(z) ~ x * y - beta * z]
sys = ODESystem(eqs; tspan=(0, 100), name=:lorenz)
prob = ODEProblem(sys, [x => 40], (0, 50), [σ => 0]) # σ is not a parameter, lets be stricter when checking that
@profview ODEProblem(sys, [x => 40], (0, 50), [sig => 0])

