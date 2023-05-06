using ExponentialNoiseSDE
using DelimitedFiles
using Test
include("small_data_test.jl")


## Small data test
small_data_test()

## Validate calculations (against MATLAB calculations)

X_data = DelimitedFiles.readdlm("./data/exampleData.txt")[:]

obs = Observation(
    X_data,
    0.0050
)

## Conditional moment settings
nTimeShiftSamplePoints = 15
nEvalPoints = 20
xEvalLims = (-1.0,1.0)
kernel = "Epanechnikov"
bandwidth = 0.1

# Variables/attributes
timeShiftSamplePoints = collect(1:nTimeShiftSamplePoints)
momentSize = (nTimeShiftSamplePoints, nEvalPoints)

momentSettings = ConditionalMomentSettings(
    timeShiftSamplePoints,
    nEvalPoints,
    xEvalLims,
    kernel,
    bandwidth,
)

# Conditional moments
conditionalMoments = build_moments(obs, momentSettings)

## Model estimates
# Settings
modelSettings = ModelEstimateSettings()

# Estimating
modelEstimate = estimate_model(conditionalMoments, modelSettings)

## Load MATLAB calculations

xEvalPointsMATLAB = DelimitedFiles.readdlm("./data/evalPointsMATLAB.txt")[:]

@show xEvalPointsMATLAB