## Calibrate tolerances for test/scientific.jl
# Runs the recovery, simulator and shrinkage checks over many seeds and prints
# tolerance = 1.5 × 99th percentile of the observed error. Rerun this when the
# estimator changes, and update the tolerances in test/scientific.jl.
#
# From the repository root:
#   julia -e 'using Pkg; Pkg.activate(temp=true); Pkg.develop(path="."); Pkg.add("StableRNGs");
#             include("test/calibration/calibrate.jl")'
# Set NSEEDS to change the number of seeds (default 100).
using ExponentialNoiseSDE, StableRNGs, Statistics
include(joinpath(@__DIR__, "..", "simulate.jl"))
include(joinpath(@__DIR__, "..", "helpers.jl"))

const θTRUE = 0.01
const DT = 0.005
const N = 50_000
const nSeeds = parse(Int, get(ENV, "NSEEDS", "100"))

f, g = example_drift, example_noise
θErr = Float64[]; fErr = Float64[]; gErr = Float64[]
nDiverged = 0
for seed in 1:nSeeds
    X = try
        simulate_process(StableRNG(seed), f, g, θTRUE, DT, N)
    catch e
        e isa DomainError || rethrow()
        global nDiverged += 1
        continue
    end
    m = fit_model(X, DT)
    x = eval_points(m); inner = abs.(x) .<= 0.6
    push!(θErr, abs(correlation_time(m) - θTRUE)/θTRUE)
    push!(fErr, maximum(abs, drift(m)[inner] .- f.(x[inner])))
    push!(gErr, maximum(abs, noise(m)[inner] .- g.(x[inner])))
end
q(v) = quantile(v, 0.99)
println("Diverged simulations skipped: $nDiverged")
println("Recovery over $(length(θErr)) seeds (N=$N, θ=$θTRUE, dt=$DT):")
println("  θ relative error : median $(round(median(θErr), sigdigits=3)), 99% $(round(q(θErr), sigdigits=3)) -> rtol $(round(1.5q(θErr), sigdigits=2))")
println("  f max abs error  : median $(round(median(fErr), sigdigits=3)), 99% $(round(q(fErr), sigdigits=3)) -> atol $(round(1.5q(fErr), sigdigits=2))")
println("  g max abs error  : median $(round(median(gErr), sigdigits=3)), 99% $(round(q(gErr), sigdigits=3)) -> atol $(round(1.5q(gErr), sigdigits=2))")
# Seed 1 is the seed used by test/scientific.jl
println("  test seed 1      : θ $(round(θErr[1], sigdigits=3)), f $(round(fErr[1], sigdigits=3)), g $(round(gErr[1], sigdigits=3))")

# Correlation time from every other lag
θErrSparse = Float64[]
for seed in 1:nSeeds
    X = try
        simulate_process(StableRNG(seed), f, g, θTRUE, DT, N)
    catch e
        e isa DomainError || rethrow()
        continue
    end
    push!(θErrSparse, abs(correlation_time(fit_model(X, DT; lags = 2:2:30)) - θTRUE)/θTRUE)
end
println("  θ relative error, lags 2:2:30: median $(round(median(θErrSparse), sigdigits=3)), 99% $(round(q(θErrSparse), sigdigits=3)) -> rtol $(round(1.5q(θErrSparse), sigdigits=2)); seed 1: $(round(θErrSparse[1], sigdigits=3))")

# Linear process with constant noise: cannot diverge, and has the analytic
# stationary variance Var X = g²/(2a(1 + aθ)) (checks the η normalisation)
const A = 2.5
fLin(x) = -A*x
gLin(x) = 1.0
varTrue = 1/(2A*(1 + A*θTRUE))
varErr = [abs(var(simulate_process(StableRNG(s), fLin, gLin, θTRUE, DT, 200_000)) - varTrue)/varTrue for s in 1:nSeeds]
println("Linear process Var X relative error (N=200000): median $(round(median(varErr), sigdigits=3)), 99% $(round(q(varErr), sigdigits=3)) -> rtol $(round(1.5q(varErr), sigdigits=2)); seed 1: $(round(varErr[1], sigdigits=3))")

# Spread of θ̂ across seeds should shrink like 1/sqrt(N)
spread(n, seeds) = std([correlation_time(fit_model(simulate_process(StableRNG(s), fLin, gLin, θTRUE, DT, n), DT)) for s in seeds])
ratios = [spread(160_000, (8k+1):(8k+8)) / spread(10_000, (8k+1):(8k+8)) for k in 0:9]
println("Spread ratio std(θ̂ | 16N)/std(θ̂ | N), 10 groups of 8 seeds: ", round.(ratios, digits=3), "  (expect ≈ 0.25)")
