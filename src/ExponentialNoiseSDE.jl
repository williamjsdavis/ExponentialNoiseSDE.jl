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

export ModelEstimateSettings, ModelEstimate
export estimate_model

include("utils.jl")

include("Observation.jl")

include("kernels.jl")
include("KBR.jl")
include("ConditionalMoments.jl")

include("autocorrelation.jl")
include("ModelEstimate.jl")

end # module ExponentialNoiseSDE
