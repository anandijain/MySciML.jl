using ModelingToolkit
using ModelingToolkit: ParentScope
using DifferentialEquations
using ModelingToolkit
using ModelingToolkit: ParentScope

function a(; name)
    @parameters t
    @variables x(t)[1:2], xx(t)[1:2]
    D = Differential(t)

    x = ParentScope.(Symbolics.scalarize(x))
    xx = Symbolics.scalarize(xx)

    eq = [D(xx[1]) ~ x[1], D(xx[2]) ~ x[2]]
    sts = [x; xx]
    any(hasmetadata.(sts, SymScope))

    asys = ODESystem(eq, t, sts, []; name)
    Debugger.@enter ODESystem(eq, t, sts, []; name)
    return asys
end
@named asys = a()
asts = states(asys)
any(hasmetadata.(asts, SymScope)) 

function main(; name)
    @parameters t
    @variables x(t)[1:2]
    D = Differential(t)

    @named s1 = a()
    any(hasmetadata.(s1.states, SymScope)) 
    # (; x) = s1
    eqs = [D(x[1]) ~ 1, D(x[2]) ~ 1]
    sys = compose(ODESystem(eqs; name), s1)
    return sys
end

@named sys = main()
structural_simplify(sys)

function a(; name)
    @parameters t
    @variables x(t)[1:2] xx(t)[1:2] = 0
    x = Symbolics.scalarize(x)
    xx = Symbolics.scalarize(xx)
    D = Differential(t)
    eq = [D(xx[1]) ~ ParentScope(x[1]), D(xx[2]) ~ ParentScope(x[2])]
    return ODESystem(eq; name)
end

function main(; name)
    @parameters t
    @variables x(t)[1:2] = 0
    D = Differential(t)


    @named s1 = a()
    eqs = [D(x[1]) ~ 1, D(x[2]) ~ 1]
    sys = compose(ODESystem(eqs; name), s1)
    return sys
end

@named sys = main()
ssys = structural_simplify(sys)
prob = ODEProblem(ssys, [], (0, 100))
sol = solve(prob)

eqs = [
    D(xx[1]) ~ D(x[1])
    D(x[1]) ~ 1
]
@named s2 = ODESystem(eqs)
ssys = structural_simplify(s2)



function b(; name)
    @parameters t
    @variables x(t) xx(t) = 0

    D = Differential(t)
    eq = [D(xx) ~ ParentScope(x)]
    return ODESystem(eq; name)

end

function mainb(; name)
    @parameters t
    @variables x(t) = 0
    D = Differential(t)


    @named s1 = b()
    eqs = [D(x) ~ 1]
    sys = compose(ODESystem(eqs; name), s1)
    return sys
end

@named sysb = mainb()
equations(sysb)