module ExponentialNoiseSDE

greet() = print("Hello World!")

export Observation

export ConditionalMoments, ConditionalMomentSettings
export Epanechnikov
export build_moments
export ModelEstimateSettings, ModelEstimate

include("Observation.jl")

include("kernels.jl")
include("KBR.jl")
include("ConditionalMoments.jl")

#include("autocorrelation.jl")
include("ModelEstimate.jl")

end # module ExponentialNoiseSDE
