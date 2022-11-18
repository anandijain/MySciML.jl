module MySciML

# using SciMLBase, DiffEqBase, OrdinaryDiffEq, DifferentialEquations, ModelingToolkit, Symbolics
using Pkg
paths = [
    "/Users/anand/.julia/dev/SciMLBase",
    "/Users/anand/.julia/dev/DiffEqBase",
    "/Users/anand/.julia/dev/OrdinaryDiffEq",
    "/Users/anand/.julia/dev/DifferentialEquations",
    "/Users/anand/.julia/dev/ModelingToolkit",
    "/Users/anand/.julia/dev/Symbolics"]
function upit(; pull=true)
    # if pull
    for path in paths
        upit(path;pull)
    end
    # end
    # Pkg.update()
end

function upit(path;pull=true)
    cd(path)
    pull && run(`git pull`)

    Pkg.activate(path)
    # Pkg.instantiate()
    Pkg.update()
    cd(@__DIR__)
    Pkg.activate(dirname(@__DIR__))
end

export upit

end # module MySciML



