using ModelingToolkit, MethodOfLines, OrdinaryDiffEq, DomainSets

@parameters x y t
@variables u(..) v(..)
Dt = Differential(t)
Dx = Differential(x)
Dy = Differential(y)
Dxx = Differential(x)^2
Dyy = Differential(y)^2

∇²(u) = Dxx(u) + Dyy(u)

brusselator_f(x, y, t) = (((x - 0.3)^2 + (y - 0.6)^2) <= 0.1^2) * (t >= 1.1) * 5.0

x_min = y_min = t_min = 0.0
x_max = y_max = 1.0
t_max = 11.5

α = 10.0

u0(x, y, t) = 22(y * (1 - y))^(3 / 2)
v0(x, y, t) = 27(x * (1 - x))^(3 / 2)

eq = [Dt(u(x, y, t)) ~ 1.0 + v(x, y, t) * u(x, y, t)^2 - 4.4 * u(x, y, t) + α * ∇²(u(x, y, t)) + brusselator_f(x, y, t),
    Dt(v(x, y, t)) ~ 3.4 * u(x, y, t) - v(x, y, t) * u(x, y, t)^2 + α * ∇²(v(x, y, t))]

domains = [x ∈ Interval(x_min, x_max),
    y ∈ Interval(y_min, y_max),
    t ∈ Interval(t_min, t_max)]

# Periodic BCs
bcs = [u(x, y, 0) ~ u0(x, y, 0),
    u(0, y, t) ~ u(1, y, t),
    u(x, 0, t) ~ u(x, 1, t), v(x, y, 0) ~ v0(x, y, 0),
    v(0, y, t) ~ v(1, y, t),
    v(x, 0, t) ~ v(x, 1, t)]

@named pdesys = PDESystem(eq, bcs, domains, [x, y, t], [u(x, y, t), v(x, y, t)])

N = 32

dx = (x_max - x_min) / N
dy = (y_max - y_min) / N

order = 2

discretization = MOLFiniteDifference([x => dx, y => dy], t, approx_order=order, grid_align=center_align)
prob = discretize(pdesys, discretization)
sol = solve(prob, TRBDF2(), saveat=0.1)