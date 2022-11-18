using CSV, DataFrames

using Symbolics, ModelingToolkit

@parameters t sig = 10 rho = 28.0 beta = 8 / 3    #= /Users/anand/.julia/config/startup.jl:260 =#
@variables x(t) = begin    #= /Users/anand/.julia/config/startup.jl:261 =#
    100
end y(t) = begin
    1.0
end z(t) = begin
    1
end
D = Differential(t)
eqs = [D(x) ~ sig * (y - x), D(y) ~ x * (rho - z) - y, D(z) ~ x * y - beta * z]
sys = ODESystem(eqs; name=:lorenz)

cnames = subtypes(Symbolics.AbstractVariableMetadata)
sps = ModelingToolkit.SYMBOLIC_METADATA_PAIRS
pushfirst!()
a, b = first.(sps), last.(sps)

last_ones = [
    :default => :VariableDefault
]
a .=> getfield.((ModelingToolkit,), last.(sps))

# best order
# everything but name and context is optional 
cnames = [
    :name, # string (fukc it we eval) 
    :context, # symbol ?
    :default, # number?
    :unit, # string
    :description, # string (could also be reference to publication)
    :bounds, # 2-tuple of numbers for now 
    :tunable, # bool 
    :dist, # distribution
    :connect, # bool
    :noise, # 
    :input,
    :output,
    :irreducible,
    :state_priority,
    :disturbance,
    :binary,
    :integer,
]

ModelingToolkit.MTKConstantCtx
ModelingToolkit.MTKParameterCtx
ModelingToolkit.VariableBinary
ModelingToolkit.VariableBounds
ModelingToolkit.VariableConnectType
ModelingToolkit.VariableDescription
ModelingToolkit.VariableDistribution
ModelingToolkit.VariableDisturbance
ModelingToolkit.VariableInput
ModelingToolkit.VariableInteger
ModelingToolkit.VariableIrreducible
ModelingToolkit.VariableNoiseType
ModelingToolkit.VariableOutput
ModelingToolkit.VariableStatePriority
ModelingToolkit.VariableTunable
ModelingToolkit.VariableUnit
Symbolics.VariableDefaultValue
Symbolics.VariableSource

df = DataFrame(cnames .=> Ref(Any[]))

using Similitude

length(Metric) / time(Metric)
