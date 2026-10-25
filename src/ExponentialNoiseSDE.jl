module ExponentialNoiseSDE

#TODO: update code to use mean()
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
