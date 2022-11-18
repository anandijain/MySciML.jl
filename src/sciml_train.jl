using Distributions, DiffEqFlux, Noise, Plots, StatsPlots, Flux, Optimization, DifferentialEquations, CSV, DataFrames
af = CSV.read("dataset1.csv", DataFrame, header = false)
data = Array{Float32}(af)           
# The dataset contrains two group of measurement. 
# The first three rows are measurement 1 and the last three rows are measurement 2.
# The first column contains intial states for two measurement 

n_exp = 200;

dudt2 = Chain(Dense(3,16,relu),
           Dense(16,23,relu),
           Dense(23,16,relu),
           Dense(16,3))           

tspan0 = range(0.0f0,0.05,5)
tspan1 = range(0.051f0,0.19,5)
tspan2 = range(0.2,2.0f0,20)
tspan = (0.0f0, 2.0f0)
tsteps = vcat(tspan0, tspan1, tspan2)
para = [2.4e4, 0.15, 0.18] 

prob_neuralode = NeuralODE(dudt2, tspan, Tsit5(), saveat = tsteps)


function loss_neuralode(p)
    loss = 0.0f0;
    for i = 1:2
        index = (i-1)*3+1;        
        ym = data[index:index+2, :];    # 30*3
        pred = Array(prob_neuralode(data[index:index+2, 1], p))   
        loss = loss + sum(abs2, ym .- pred)   # + 2 * sum(abs2, ym[:, 2] .- nny[:, 2])   # I gave
        return loss
    end
  return loss #, pred2 , pred3 , pred4 
end

callback = function (p, l, pred; doplot = true)
    display(l)
    return false
end
  
# using Optimization
result_neuralode = DiffEqFlux.sciml_train(loss_neuralode, prob_neuralode.p,
                                            ADAM(0.05), cb = callback,
                                            maxiters = 2000)