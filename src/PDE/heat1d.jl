using DifferentialEquations, ModelingToolkit, MethodOfLines, DomainSets
using Plots

# Method of Manufactured Solutions: exact solution

# ModelingToolkit.parameters(sys::PDESystem) = getfield(sys, :ps)
# ModelingToolkit.states(sys::PDESystem) = getfield(sys, :dvs)
# ModelingToolkit.get_states(sys::PDESystem) = getfield(sys, :dvs)
# ModelingToolkit.has_states(sys::PDESystem) = isdefined(sys, :dvs)
# ModelingToolkit.has_ps(sys::PDESystem) = isdefined(sys, :ps)

u_exact = (x, t) -> exp.(-t) * cos.(x)
function heat1d(xmin, xmax, tmin, tmax; name)
    @parameters t x
    @variables u(..)
    Dt = Differential(t)
    Dxx = Differential(x)^2

    eq = Dt(u(t, x)) ~ Dxx(u(t, x))

    bcs = [u(tmin, x) ~ cos(x),
        u(t, xmin) ~ exp(-t),
        u(t, xmax) ~ exp(-t) * cos(xmax)]

    domains = [t ∈ Interval(tmin, tmax),
        x ∈ Interval(xmin, xmax)]
    PDESystem(eq, bcs, domains, [t, x], [u(t, x)], []; name)
end

@named sys1 = heat1d(0, 1, 0, 10)
@named sys1 = heat1d(1, 2, 0, 10)

@variables x t
@variables u(..)
@parameters a b
Dx = Differential(x)

# connect_eqs = [
#     getproperty(sys1, :u)(t, 1) ~ sys2.u(t, 1),
#     a * Dx(x(t, 1)) - b * Dx(x(t, 2)) ~ 0
# ]

dx = 0.1
order = 2
discretization = MOLFiniteDifference([x => dx], t)
@named pdesys = heat1d(0, 1, 0, 10)
prob = discretize(pdesys, discretization)
sol = solve(prob, Tsit5(), saveat=0.2)

discrete_x = sol[x]
discrete_t = sol[t]
solu = sol[u(t, x)]

plt = plot()

for i in 1:length(discrete_t)
    plot!(discrete_x, solu[i, :], label="Numerical, t=$(discrete_t[i])")
    scatter!(discrete_x, u_exact(discrete_x, discrete_t[i]), label="Exact, t=$(discrete_t[i])")
end
display(plt)
savefig("plot.png")


@register_symbolic foo(t, x)
@register_symbolic foo(t=1, x=2)




# CallWithIvs

@variables t u(; t, x) #u(;t, x[1:5])
@variables x[1:5] u(t, x) u(t, x)[1:3]


u => u(t, x)
u(t=1) => u(1, x)

@parameters t x[1:5]
@variables u(..)

# ??
Differential(t)(x)
collect(Differential(t).(x))


@variables t x[1:3] u(t, x)
tdom = t ∈ Interval(0 .. 1)
xdoms = [xi ∈ Interval(0 .. 1) for xi in x]
doms = [tdom; xdoms]

xdoms2 = x ∈ Ball()


ndim(Sphere())
# (Disk()
eval(lorenz())
prob = ODEProblem(sys)
sol = solve(prob)
sol(0) == sol[1]
@test_throws sol(0, 1)

abstract type CallW end

struct CallWithKwargs
    f
    args
end

(x::CallWithKwargs)(; args) = x.f(x.args..., args...)