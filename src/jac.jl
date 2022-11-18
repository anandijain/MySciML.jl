using SciMLBase, DiffEqBase, OrdinaryDiffEq, DifferentialEquations, ModelingToolkit, Symbolics
using Cthulhu
function lorenz()
    :(
        begin
            @parameters t sig = 10 rho = 28.0 beta = 8 / 3
            @variables x(t) = 100 y(t) = 1.0 z(t) = 1
            D = Differential(t)

            eqs = [D(x) ~ sig * (y - x),
                D(y) ~ x * (rho - z) - y,
                D(z) ~ x * y - beta * z]
            sys = ODESystem(eqs; name=:lorenz)  
        end
    )
end

eval(lorenz())
sys2 = deepcopy(sys)

@test sys2 !== sys

prob = ODEProblem(sys, [], (0, 10))

ModelingToolkit.calculate_jacobian(sys2)
prob2 = ODEProblem(sys2, [], (0, 10))
prob3 = ODEProblem(sys2, [], (0, 10);jac=true)

f = ODEFunction(sys2)
@test SciMLBase.__has_jac(f)
# @descend ODEFunction(sys2; jac=true)

@test isnothing(prob2.f.jac)
f = ODEFunction(sys2; jac=true)
alg = FBDF()