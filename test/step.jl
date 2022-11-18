using ModelingToolkit, DifferentialEquations
SR = 48000
dt = 1 // SR
@parameters t
@variables x(t) = 1.0 vx = 0.0
D = Differential(t)
eq = [D(D(x)) ~ -1000^2 * x]
@named harmosc = ODESystem(eq)
ssys = structural_simplify(harmosc)
prob = ODEProblem(ssys, [D(x) => 0], (0, 1000); dt=dt, save_everystep=false)
integ = init(prob, Tsit5(); abstol=1e-6, reltol=1e-6, dt=dt, save_everystep=false, verbose=false)
integ = init(prob, Tsit5(); abstol=1e-6, reltol=1e-6, dt=dt, save_everystep=false)
integ = init(prob, Tsit5(); abstol=1e-6, reltol=1e-6, dt=dt, save_everystep=false, force_dtmin=true)

integ = init(prob, Tsit5(); abstol=1e-6, reltol=1e-6, dt=dt, save_everystep=false, force_dtmin=true, maxiters=typemax(Int))
integ.do_error_check = false
solve!(integ)

for i in 1:integ.opts.maxiters+1
    step!(integ, dt, true)
end

tt = integ.t
step!(integ, dt, true)
@test integ.t == tt
integ.do_error_check = false
step!(integ, dt, true)
@test_broken integ.t == tt
@test integ.do_error_check

integ.t
for i in take(integ, buf_size)
    buf[i] = clamp(integ.u[1], -1, 1)
end
