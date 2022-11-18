using ModelingToolkit#, DifferentialEquations, Plots
# using TimerOutputs, StatsBase
# using WAV

@info "usings"

function Kuramoto(; N, natural_frequencies=ones(N), thetas=0, coupling_strength=0.7, name)
    @parameters coupling_strength = coupling_strength
    @parameters natural_frequencies[1:N] = natural_frequencies
    @variables theta(t)[1:N] = thetas

    eqs = []
    for i in 1:N
        term = coupling_strength / N * sum(sin(theta[j] - theta[i]) for j in 1:N)
        eq = D(theta[i]) ~ natural_frequencies[i] + term
        push!(eqs, eq)
    end
    ODESystem(eqs; name=name)
end


function Kuramoto2(; N, natural_frequencies=ones(N), thetas=0, coupling_strength=0.7, name)
    @parameters coupling_strength = coupling_strength
    @parameters natural_frequencies[1:N] = natural_frequencies
    @variables theta(t)[1:N] = thetas

    eqs = []
    for i in 1:N
        ti = theta[i] % 2pi
        s = 0
        for j in 1:N
            tj = theta[j] % 2pi
            s += sin(tj - ti)
        end
        term = s * coupling_strength / N
        eq = D(theta[i]) ~ natural_frequencies[i] + term
        push!(eqs, eq)
    end
    ODESystem(eqs; name=name)
end

"kura3 has unique coupling strength for each ij pair (asymetric) "
function Kuramoto3(; N, natural_frequencies=ones(N), thetas=0, coupling_strengths=0.7, name)
    @parameters coupling_strengths[1:N, 1:N] = coupling_strengths
    @parameters natural_frequencies[1:N] = natural_frequencies
    @variables theta(t)[1:N] = thetas

    eqs = []
    for i in 1:N
        ti = theta[i] % 2pi
        s = 0
        for j in 1:N
            tj = theta[j] % 2pi
            s += coupling_strengths[i, j] * sin(tj - ti)
        end
        term = s * 1 / N
        eq = D(theta[i]) ~ natural_frequencies[i] + term
        push!(eqs, eq)
    end
    ODESystem(eqs; name=name)
end
"sde"
function Kuramoto4(; N, tspan, noise_coef=0.1, natural_frequencies=ones(N), thetas=0, coupling_strengths=0.7, name)
    @parameters coupling_strengths[1:N, 1:N] = coupling_strengths
    @parameters natural_frequencies[1:N] = natural_frequencies
    @variables theta(t)[1:N] = thetas

    eqs = []
    neqs = []
    for i in 1:N
        ti = theta[i] % 2pi
        s = 0
        for j in 1:N
            tj = theta[j] % 2pi
            s += coupling_strengths[i, j] * sin(tj - ti)
        end
        term = s * 1 / N
        eq = D(theta[i]) ~ natural_frequencies[i] + term
        push!(neqs, theta[i] * randn() * noise_coef)
        push!(eqs, eq)
    end
    SDESystem(ODESystem(eqs; name), neqs; tspan, name)
end

"stochastic on n-dims with"
function Kuramoto5(; N, tspan, dim::Integer, noise_coef=0.1, natural_frequencies=ones(N), thetas=0, coupling_strengths=0.7, name)
    @parameters coupling_strengths[1:N, 1:N, 1:dim] = coupling_strengths
    @parameters natural_frequencies[1:N, 1:dim] = natural_frequencies
    @variables theta(t)[1:N, 1:dim] = thetas

    eqs = []
    neqs = []
    for k in 1:dim
        for i in 1:N
            ti = theta[i, k] % 2pi
            s = 0
            for j in 1:N
                tj = theta[j, k] % 2pi
                s += coupling_strengths[i, j, k] * sin(tj - ti)
            end
            term = s * 1 / N
            eq = D(theta[i, k]) ~ natural_frequencies[i, k] + term
            push!(eqs, eq)
            push!(neqs, theta[i, k] * randn() * noise_coef)
        end
    end
    # end
    sys = ODESystem(eqs; tspan, name)
    SDESystem(sys, neqs; name)
end

"n-dims"
function Kuramoto6(; N, tspan, dim::Integer, noise_coef=0.1, natural_frequencies=ones(N), thetas=0, coupling_strengths=0.7, name)
    @parameters coupling_strengths[1:N, 1:N, 1:dim] = coupling_strengths
    @parameters natural_frequencies[1:N, 1:dim] = natural_frequencies
    @variables theta(t)[1:N, 1:dim] = thetas

    eqs = []
    neqs = []
    for k in 1:dim
        for i in 1:N
            ti = theta[i, k] % 2pi
            s = 0
            for j in 1:N
                tj = theta[j, k] % 2pi
                s += coupling_strengths[i, j, k] * sin(tj - ti)
            end
            term = s * 1 / N
            eq = D(theta[i, k]) ~ natural_frequencies[i, k] + term
            push!(eqs, eq)

            neq = theta[i, k] * t/100
            push!(neqs, neq)
        end
    end
    # end
    sys = ODESystem(eqs; tspan, name)
    SDESystem(sys, neqs; name)
end


# @timeit to "$N" kuramoto = ODESystem(eqs; name=Symbol(:kuramoto, N))
# @timeit toss "$N" sys = structural_simplify(kuramoto)

# @btime ODESystem(eqs; name=Symbol(:kuramoto, N))
# @profview ODESystem(eqs; name=Symbol(:kuramoto, N))

# if the 
@parameters t
D = Differential(t)

N = 10
dim = 1
const to = TimerOutput()
const toss = TimerOutput()

u0s = 2pi * rand(N)
# kuramoto = Kuramoto(N=N, thetas=u0s, name=Symbol(:kuramoto, N))
# kuramoto = Kuramoto2(N=N, thetas=u0s, name=Symbol(:kuramoto, N))
t1 = 100
tspan = (0, t1)

coupling_strengths = randn(N, N)
natural_frequencies = rand(N, dim) * 440
name = Symbol(:kuramoto, N, dim)
thetas = u0s
# kuramoto = Kuramoto4(N=N, noise_coef=1e-2, tspan=tspan, thetas=u0s, natural_frequencies=rand(N) * 440, name=Symbol(:kuramoto, N))
kuramoto = Kuramoto6(;N, tspan, dim, natural_frequencies, thetas, name)
@info "sys"
u0 = states(kuramoto) .=> u0s
# sys = structural_simplify(kuramoto)


u0 = states(kuramoto) .=> reduce(vcat, u0s)

# sys = structural_simplify(kuramoto)
@info "ssys"
# @profview structural_simplify(kuramoto)
# end
# @named kura[1:1000] = Kuramoto(N=10)
# probs = map(x->ODEProblem(x, [], tspan), kura)
# ep = EnsembleProblem(probs)
# sol = solve.(probs)
# @which EnsembleProblem(probs)
# sol = solve(probs)

# function EnsembleProblem(Array{ODEProblem,1})

#     ODEProblem(f, u0, tspan)
# end
sys = kuramoto
# prob = ODEProblem(sys, [], tspan)
prob = SDEProblem(sys)
@info "prob"
# VSCodeServer.@enter SDEProblem(sys)
@profview SDEProblem(sys)