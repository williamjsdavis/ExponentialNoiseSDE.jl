module ExponentialNoiseSDE

greet() = print("Hello World!")

#TODO: update code to use mean()
using FFTW: fft, ifft
using NLopt: Opt, optimize
using Statistics: mean, var

export Observation

export Epanechnikov
export ConditionalMoments, ConditionalMomentSettings
export build_moments

export FiniteDiff, SmoothedFiniteDiff
export ModelEstimateSettings, ModelEstimate
export estimate_model

include("utils.jl")

include("Observation.jl")

include("kernels.jl")
include("KBR.jl")
include("ConditionalMoments.jl")

include("autocorrelation.jl")

include("ModelEstimate.jl")
include("theta_estimation.jl")
include("lambda_estimation.jl")
include("function_estimation.jl")
include("model_error.jl")
include("default_settings.jl")

end # module ExponentialNoiseSDE
