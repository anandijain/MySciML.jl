using ModelingToolkit
using Symbolics: unwrap
using SymbolicUtils
using SymbolicUtils.Code
using ModelingToolkit: isvariable
import Symbolics: value
using Symbolics, Test, ModelingToolkit
using ModelingToolkit: is_delay_var, isvariable
using DelayDiffEq
using Plots

function my_is_delay_var(iv::SymbolicUtils.Symbolic, v)
    isvariable(v) || return false
    istree(v) || return false
    if operation(v) === getindex
        v = arguments(v)[1]
    end
    istree(v) || return false
    any(!isequal(iv, arg) && occursin(iv, arg) for arg in arguments(v))
end

my_is_delay_var(ivs, v) = any((my_is_delay_var(iv, v) for iv in ivs))

function collect_delay_variables(sys)
    ivs = independent_variables(sys)
    eqs = equations(sys)
    varss = Set()
    for eq in eqs
        ModelingToolkit.vars!(varss, eq)
    end
    collect(filter(x -> my_is_delay_var(ivs, x), varss))
end

has_delay(sys) = !isempty(collect_delay_variables(sys))

@variables t x(..)[1:3]
pars = @parameters p0, q0, v0, d0, p1, q1, v1, d1, d2, beta0, beta1, tau
D = Differential(t)
hist3 = x(t - tau)[3]
x = collect(x(t))
eqs = [D(x[1]) ~ (v0 / (1 + beta0 * (hist3^2))) * (p0 - q0) * x[1] - d0 * x[1]
    D(x[2]) ~ (v0 / (1 + beta0 * (hist3^2))) * (1 - p0 + q0) * x[1] +
              (v1 / (1 + beta1 * (hist3^2))) * (p1 - q1) * x[2] - d1 * x[2]
    D(x[3]) ~ (v1 / (1 + beta1 * (hist3^2))) * (1 - p1 + q1) * x[2] - d2 * x[3]]
@named sys = ODESystem(eqs, t, x, pars)

delay_terms = collect_delay_variables(sys)
ivs = independent_variables(sys)
iv = unwrap(only(ivs))
eqs = equations(sys)
hh = Sym{Any}(:h)

# HACK
sub_term = term(getindex, term(hh, pars, unwrap(iv - tau), type=Real), 3, type=Real)

eqs = substitute.(eqs,
    (Dict(
        delay_terms[1] => sub_term),))
out = Sym{Any}(:out)
body = SetArray(false, out, getfield.(eqs, :rhs))
func = Func([out, DestructuredArgs(states(sys)), hh, DestructuredArgs(parameters(sys)), iv],
    [], body)

my_func_expr = toexpr(func)
my_func = eval(my_func_expr)
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

prob = DDEProblem(my_func, u0, h, tspan, p; constant_lags=lags)
alg = MethodOfSteps(Rodas4())
sol = solve(prob, alg)

dt = 0.1
dt2 = 0.2
@variables t x(t) y(t) u(t) r(t) yd1(t) ud1(t) yd2(t) ud2(t)
@parameters kp
D = Differential(t)

dt = 0.1
@variables t x(t) y(..) u(t) yd(t) ud(t) r(t)
@parameters kp
D = Differential(t)

eqs = [yd ~ Sample(t, dt)(y(t - 1))
    ud ~ kp * (r - yd)
    # plant (time continuous part)
    u ~ Hold(ud)
    D(x) ~ -x + u
    y(t) ~ x]

@named sys = ODESystem(eqs)
@test has_delay(sys) == true




# sol = NDSolve[{x''[t] + x[t-1] == 0, x[t/t<=0] == t^2}, x, {t, -1, 5}]


@variables x(..)
D = Differential(t)
D2 = D^2
deq = D2(x(t)) ~ -x(t - 1)
@named sys = ODESystem(deq; tspan)
ode_order_lowering(sys)
ineq = Inequality(t, 0, <=)
x(t), t < 0

@parameters t x
@variables u(..)
Dt = Differential(t)
Dx = Differential(x)
Dxx = Dx^2









function f(ddu, du, u, p, t)
    ddu[1] = -u[1]
end
du0 = [0]
u0 = [1]
harmosc = SecondOrderODEProblem(f, du0, u0, tspan)
sys = modelingtoolkitize(harmosc)
ssys = structural_simplify(sys)


using ModelingToolkit, DifferentialEquations

@parameters t
@variables x(t) = 1.0 vx = 0.0
D = Differential(t)
eq = [D(D(x)) ~ -x]
@named harmosc2 = ODESystem(eq)
sys2 = ode_order_lowering(harmosc2)
ssys2 = structural_simplify(harmosc2)
equations(ssys)
prob = ODEProblem(ssys, [D(x) => 0], (0, 10))
sol = solve(prob)
plot(sol)
@test isapprox(cos.(sol.t), sol[x]; rtol=1e-4)

@test ModelingToolkit.isisomorphic(sys, sys2)
