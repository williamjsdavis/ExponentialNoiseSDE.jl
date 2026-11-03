## Benchmark: time of each stage of the estimation against the series length
# From the repository root:
#   julia --project=benchmark -e 'using Pkg; Pkg.develop(path="."); Pkg.instantiate()'
#   julia --project=benchmark benchmark/scaling.jl
# Prints a table of times and writes benchmarking_time.svg in benchmark/.
# Only the public API is timed, so the stages are the ones a user can call.
using ExponentialNoiseSDE
using BenchmarkTools
using Plots
using StableRNGs

# Ground-truth simulator of the model, shared with the tests
include(joinpath(@__DIR__, "..", "test", "simulate.jl"))

# Series lengths, and settings of the example in example/
const LENGTHS = [50_000, 100_000, 1_000_000, 10_000_000]
const DT = 0.005
momentSettings = ConditionalMomentSettings(
    timeShiftSamplePoints = 1:15,
    nEvalPoints = 20,
    evalLims = (-1.0, 1.0),
    bandwidth = 0.1,
)
fitSettings = ModelEstimateSettings()

# Linear drift and constant noise, so long paths cannot diverge; the cost
# depends on the series length, not on the model
simulate(N) = simulate_process(StableRNG(1), x -> -2.5*x, x -> 1.0, 0.01, DT, N)

# Minimum time in seconds of each stage over 5 runs (one evaluation each, and a
# time budget large enough for the longest series); the moments include the
# autocorrelation
function stage_times(N)
    observation = Observation(simulate(N), DT)
    conditionalMoments = build_moments(observation, momentSettings)
    moments = @belapsed build_moments($observation, $momentSettings) samples=5 evals=1 seconds=120
    fit = @belapsed estimate_model($conditionalMoments, $fitSettings) samples=5 evals=1 seconds=120
    total = @belapsed estimate_model($observation, $momentSettings, $fitSettings) samples=5 evals=1 seconds=120
    return (; moments, fit, total)
end

println("ExponentialNoiseSDE ", pkgversion(ExponentialNoiseSDE), ", Julia ", VERSION,
    ", ", Sys.cpu_info()[1].model, ", ", Threads.nthreads(), " thread(s)")
times = map(stage_times, LENGTHS)

println("\n| Samples | Moments (s) | Fit from moments (s) | Total (s) |")
println("|---|---|---|---|")
for (N, t) in zip(LENGTHS, times)
    println("| $N | ", round(t.moments, sigdigits=3), " | ", round(t.fit, sigdigits=3),
        " | ", round(t.total, sigdigits=3), " |")
end

p = plot(xscale=:log10, yscale=:log10, xlabel="Number of time-series samples",
    ylabel="Time used, seconds", yticks=10.0 .^ (-4:1), ylims=(1E-4, 10),
    legend=:topleft, size=(700, 450), margin=5Plots.mm)
for (stage, label) in ((:total, "Total"),
        (:moments, "Conditional moments and autocorrelation"),
        (:fit, "Fit from moments: θ, basis functions, drift and noise"))
    plot!(p, LENGTHS, getproperty.(times, stage), marker=:circle, markersize=6, label=label)
end
savefig(p, joinpath(@__DIR__, "benchmarking_time.svg"))
