"""
    ExponentialNoiseSDE

Estimate a stochastic differential equation driven by exponentially correlated
(Ornstein-Uhlenbeck) noise from a sampled time series,

    dX/dt = f(X) + g(X) η(t),    dη = -η/θ dt + (1/θ) dW,

where f is the drift, g the noise amplitude and θ the correlation time of the
noise. η is normalised so that <η(t)η(s)> = exp(-|t-s|/θ)/(2θ): as θ → 0 it
becomes white noise and the model becomes dX = f dt + g ∘ dW (Stratonovich).

The method (Lehle & Peinke, 2018) estimates conditional moments of the
increments by kernel regression ([`build_moments`](@ref)), fits θ to the
autocorrelation increments, and solves for f and g ([`estimate_model`](@ref)).
Uncertainties come from a block bootstrap ([`estimate_bootstrap_statistics`](@ref)).
"""
module ExponentialNoiseSDE

using FFTW: fft, ifft
using NLopt: Opt, optimize
using Statistics: mean, var, std, median, quantile
using Random: rand, AbstractRNG, default_rng
using Logging: Logging, AbstractLogger, current_logger, with_logger

export Observation

export Epanechnikov
export ConditionalMoments, ConditionalMomentSettings
export build_moments

export EstimateTheta, FixedTheta
export FiniteDiff, SmoothedFiniteDiff
export ModelEstimateSettings, ModelEstimate
export estimate_model

export BootstrapSettings
export estimate_bootstrap_statistics

export correlation_time, drift, noise, eval_points

include("utils.jl")

include("Observation.jl")

include("kernels.jl")
include("KBR.jl")
include("ConditionalMoments.jl")

include("autocorrelation.jl")

include("function_estimation.jl")
include("theta_estimation.jl")
include("lambda_estimation.jl")
include("model_error.jl")
include("ModelEstimate.jl")

include("Bootstrap.jl")
include("BootstrapStatistics.jl")

end # module ExponentialNoiseSDE
