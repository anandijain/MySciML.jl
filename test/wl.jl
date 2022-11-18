using ModelingToolkit, DifferentialEquations
using ModelingToolkit
using Symbolics: unwrap
using SymbolicUtils
using SymbolicUtils.Code
using ModelingToolkit: isvariable
import Symbolics: value
using Symbolics, Test, ModelingToolkit
using ModelingToolkit: is_delay_var, isvariable
@parameters t p
@variables x(..)

eqs = [
    D(x(t)) ~ -p * x(t),
    x(0) ~ 1.0,
    p ~ 1
]

# function my_odesys(eqs)

# end
varss = Set()
for eq in eqs
    ModelingToolkit.vars!(varss, eq)
end
vs = collect(value.(varss))

for v in vs
    
    