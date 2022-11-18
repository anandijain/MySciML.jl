using ModelingToolkit, DifferentialEquations, Symbolics# , Plots
using TimerOutputs
@parameters t
D = Differential(t)

@parameters coupling_strength = 0.7
N = 50
const to = TimerOutput()
const toss = TimerOutput()
# for N in 10:50:110
# @info N
# @profview begin 
function Kuramoto(; N, coupling_strength=0.7, name)
    @parameters natural_frequencies[1:N] = rand(N)
    @variables theta(t)[1:N] = 2pi * rand(N)

    eqs = []
    for i in 1:N
        term = coupling_strength / N * sum(sin(theta[j] - theta[i]) for j in 1:N)
        push!(eqs, D(theta[i]) ~ natural_frequencies[i] + term)
    end
    ODESystem(eqs; name=name)
end
# @timeit to "$N" kuramoto = ODESystem(eqs; name=Symbol(:kuramoto, N))
# @timeit toss "$N" sys = structural_simplify(kuramoto)

# @btime ODESystem(eqs; name=Symbol(:kuramoto, N))
# @profview ODESystem(eqs; name=Symbol(:kuramoto, N))
kuramoto = Kuramoto(N=N, name=Symbol(:kuramoto, N))
sys = structural_simplify(kuramoto)
@profview structural_simplify(kuramoto)
# end
@named kura[1:1000] = Kuramoto(N=10)
t1 = 100
tspan = (0, t1)
# function EnsembleProblem()
prob = ODEProblem(sys, [], tspan)

n_steps = 1000
ts = range(0, t1, length=n_steps)
sol = solve(prob)
xs = sol(ts)
ys = xs .% 2pi # VectorOfArray

arr = Array(xs)

sts = states(sys)
anim = @animate for i in 1:size(arr, 2)
    ti = ts[i]
    plt = plot(; ticks=false, legend=false, axis=false, title="$(round(Int, ti))", xlims=(-1.2, 1.2), ylims=(-1.2, 1.2), aspect_ratio=:equal)
    last_ten = maximum([1, i - 10]):i
    nlast = length(last_ten)
    for st in sts
        thetas = xs[st][last_ten]
        a = sincos.(thetas)
        x, y = first.(a), last.(a)
        plot!(plt, x, y)
        scatter!(plt, x[end:end], y[end:end])
    end
end
mp4(anim, "kuramoto.mp4", fps=60)
# mp4 = mp4("kuramoto.mp4", anim)
@btime plot(; legend=false, axis=false, title="$(round(Int, 4.3))", xlims=(-1.2, 1.2), ylims=(-1.2, 1.2), aspect_ratio=:equal)
plot(sol)

df = timerout_to_df(to)
df2 = timerout_to_df(toss)

sort!(df, :name)
sort!(df2, :name)

plt = plot(df.name, df.time ./ 1e9)
plt = plot(df2.name, df2.time ./ 1e9)
plot!(plt, df2.name, df2.time)

@btime plot(df.name, df.time)


