
using DifferentialEquations, ModelingToolkit, MethodOfLines, DomainSets

N = 3
xmin = 0
xmax = 1
tmin = 0
tmax = 1

@parameters t x[1:N]
@variables u(..)
Dt = Differential(t)
domains = [t ∈ Interval(tmin, tmax)]
append!(domains, collect(x) .∈ Interval(xmin, xmax))

dxs = Differential.(collect(x))
utx = u(t, x)
lap = sum(map(x -> x(utx), dxs .^ 2))
eq = Dt(utx) ~ lap


# constant_bcs
bcs = [
    u(tmin, x[1], x[2]) ~ cos(x[1]) + cos(x[2]),
    u(t, xmin, x[2]) ~ cos(xmin) + cos(x[2]),
    u(t, xmax, x[2]) ~ cos(xmax) + cos(x[2]),
    u(t, x[1], xmin) ~ cos(x[1]) + cos(xmin),
    u(t, x[1], xmax) ~ cos(x[1]) + cos(xmax),
]

@named sys = PDESystem(eq, bcs, domains, [t, x[1], x[2]], [u(t, x[1], x[2])], [])

dx = 0.1
order = 2
discretization = MOLFiniteDifference(x .=> dx, t)
prob = discretize(sys, discretization)
sol = solve(prob);




# @test_throws MethodError x .∈ Interval(xmin, xmax)
# xdoms = collect(x) .∈ Interval(xmin, xmax)
# what should be done here 
# unsat_bcs = [
#     u(0, x) ~ 1,
#     u(t, 0) ~ 0
# ]





f = cos(x[1]) + cos(x[2])
fx(x) = cos(x[1]) + cos(x[2])
ui = 0:0.1:1
points = collect(Iterators.product(ui, ui))
heatmap(map(fx, points))