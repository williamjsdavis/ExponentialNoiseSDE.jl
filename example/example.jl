## Example: estimate the model for the example data, with bootstrap uncertainties
# From the repository root:
#   julia --project=example -e 'using Pkg; Pkg.develop(path="."); Pkg.instantiate()'
#   julia --project=example example/example.jl
# Writes the figures timeseries.svg, estimates.svg and bootstrap.svg in example/.
using ExponentialNoiseSDE
using DelimitedFiles
using Plots
using Plots.Measures
using Random

# Data simulated with dt = 0.005, θ = 0.01 and these drift and noise functions
X = vec(readdlm(joinpath(@__DIR__, "..", "test", "data", "exampleData.txt")))
observation = Observation(X, 0.005)
drift_true(x) = -2.5*x
noise_true(x) = x > 0 ? 1 + x^2/4 : 1.0

p = plot((1:10000) .* observation.dt, X[1:10000], label="", xlabel="t", ylabel="X(t)")
savefig(p, joinpath(@__DIR__, "timeseries.svg"))

## Model estimate
momentSettings = ConditionalMomentSettings(
    timeShiftSamplePoints = 1:15,
    nEvalPoints = 20,
    evalLims = (-1.0, 1.0),
    bandwidth = 0.1,
)
modelEstimate = estimate_model(observation, momentSettings, ModelEstimateSettings())
display(modelEstimate)
println()

x = eval_points(modelEstimate)
p1 = scatter(x, drift(modelEstimate), label="Estimate", title="Drift function", xlabel="x", ylabel="f(x)")
plot!(p1, drift_true, extrema(x)..., label="True")
p2 = scatter(x, noise(modelEstimate), label="", title="Noise function", xlabel="x", ylabel="g(x)")
plot!(p2, noise_true, extrema(x)..., label="", ylim=(0.75, Inf))
savefig(plot(p1, p2, size=(900, 400), bottom_margin=5mm, left_margin=5mm),
    joinpath(@__DIR__, "estimates.svg"))

## Bootstrap uncertainties
bootstrapStatistics = estimate_bootstrap_statistics(
    modelEstimate, BootstrapSettings(nSamples=200); rng=Xoshiro(1))
println(bootstrapStatistics)

function plot_with_interval(summary, estimate, truth; legend=true, kwargs...)
    p = plot(x, summary.lower95, fillrange=summary.upper95, linealpha=0, fillalpha=0.2,
        color=5, label=legend ? "95% interval" : "")
    scatter!(p, x, estimate, color=1, label=legend ? "Best estimate" : "")
    plot!(p, truth, extrema(x)..., color=2, label=legend ? "True" : ""; kwargs...)
    return p
end
p1 = plot_with_interval(drift(bootstrapStatistics), drift(modelEstimate), drift_true;
    title="Drift function", xlabel="x", ylabel="f(x)")
p2 = plot_with_interval(noise(bootstrapStatistics), noise(modelEstimate), noise_true;
    legend=false, title="Noise function", xlabel="x", ylabel="g(x)")
savefig(plot(p1, p2, size=(900, 400), bottom_margin=5mm, left_margin=5mm),
    joinpath(@__DIR__, "bootstrap.svg"))
