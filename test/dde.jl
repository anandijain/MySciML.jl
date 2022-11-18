using ModelingToolkit, DifferentialEquations

@parameters t p = 1 p2 = 1
@variables x(..) = 10.0
D = Differential(t)

eqs = [D(x(t)) ~ x * x(t)]

@named de = ODESystem(eqs; tspan=(0, 1000.0))
prob = ODEProblem(de, [x => 0], (0, 50), [p => 0])
@test prob.u0 == [10.0]

prob2 = ODEProblem(de, [x(t) => 0], (0, 50), [p2 => 0])
@test prob2.u0 == [0]

prob2 = ODEProblem(de, [x(t) => 0], (0, 50), [p => 0])


@variables x(t) = 10.0
D = Differential(t)
eqs = [D(x) ~ -p * x]
myvar = 1
@named de = ODESystem(eqs; tspan=(0, 1000.0))
prob = ODEProblem(de, [myvar => 0], (0, 50), [p => 0])
@test prob.u0 == [0]

prob = ODEProblem(de; defaults=[x => 40])
sol = solve(prob)

ModelingToolkit.defaults(prob::SciMLBase.AbstractODEProblem) = hasproperty(prob.f, :sys) ? ModelingToolkit.defaults(prob.f.sys) : error()
@parameters t sig = 10 rho = 28.0 beta = 8 / 3
@variables x(t) = 100 y(t) = 1.0 z(t) = 1
D = Differential(t)

eqs = [D(x) ~ sig * (y - x),
    D(y) ~ x * (rho - z) - y,
    D(z) ~ x * y - beta * z]
sys = ODESystem(eqs; tspan=(0, 100), name=:lorenz)
prob = ODEProblem(de, [x => 40], (0, 50), [σ => 0]) # σ is not a parameter, lets be stricter when checking that
prob = ODEProblem(de; u0map=[x => 40])

@profview ODEProblem(de, [x => 40], (0, 50), [sig => 0])

# cases 
x(t - p) # p parameter that could be affected by events
u[1](t - u[2]) # delayed by a state

@variables u(..)[1:10]



function bc_model(du, u, h, p, t)
    p0, q0, v0, d0, p1, q1, v1, d1, d2, beta0, beta1, tau = p
    hist3 = h(p, t - tau)[3]
    du[1] = (v0 / (1 + beta0 * (hist3^2))) * (p0 - q0) * u[1] - d0 * u[1]
    du[2] = (v0 / (1 + beta0 * (hist3^2))) * (1 - p0 + q0) * u[1] +
            (v1 / (1 + beta1 * (hist3^2))) * (p1 - q1) * u[2] - d1 * u[2]
    du[3] = (v1 / (1 + beta1 * (hist3^2))) * (1 - p1 + q1) * u[2] - d2 * u[3]
end

h(p, t) = ones(3)

tau = 1
lags = [tau]

p0 = 0.2;
q0 = 0.3;
v0 = 1;
d0 = 5;
p1 = 0.2;
q1 = 0.3;
v1 = 1;
d1 = 1;
d2 = 1;
beta0 = 1;
beta1 = 1;
p = (p0, q0, v0, d0, p1, q1, v1, d1, d2, beta0, beta1, tau)
tspan = (0.0, 10.0)
u0 = [1.0, 1.0, 1.0]

prob = DDEProblem(bc_model, u0, h, tspan, p; constant_lags=lags)
modelingtoolkitize(prob)
sol = solve(prob)
sol.alg



@variables t1
@parameters t iv2 tau
@variables u(..)[1:10]

@variables u(..)[1:10] arr(t)[1:10]
@variables x(..)
@test u(1) isa Symbolics.Term
@test u(t) isa Symbolics.Arr
@test_throws Any collect(u)
@test_throws Any length(u)
@test_throws Any length(u)

iv = value(t)

v = value(x(t - tau))
@test_broken is_delay_var(iv, v) # wrong because tau


v2 = value(x(t - 1))
@test is_delay_var(iv, v2) # correct

v3 = value(u(t - 1)[1])
is_delay_var(iv, v3)

v4 = value(u(t - 1))
@test is_delay_var(iv, v4)

v4 = value(x(iv2, t - 1))
@test_broken is_delay_var(iv, v4) # because multiple ivs

v5 = value(u(t1 - tau)[1])

v6 = value(u(iv2, t - 1))
v6 = value(u(iv2, t - 1)[1])

v7 = value(u(t - 1))


ivs = [t, iv2]
@which my_is_delay_var(iv, v)
my_is_delay_var(iv, v2)
my_is_delay_var(iv, v3)

my_is_delay_var(iv, v4)
!my_is_delay_var(iv, v5)

my_is_delay_var(iv, v6)
my_is_delay_var(iv, v6)

my_is_delay_var(iv, v7)