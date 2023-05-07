using ExponentialNoiseSDE
using DelimitedFiles
using Test
include("small_data_test.jl")

## Utilities
@testset "Statistics and utility functions" begin
    @test all(ExponentialNoiseSDE.movmean([1,2,3,4,5],3) .== [1.5,2.0,3.0,4.0,4.5])
    @test all(ExponentialNoiseSDE.movmean([1,2,3,4,5],2) .== [1.0,1.5,2.5,3.5,4.5])
    @test all(ExponentialNoiseSDE.movmean([1,2,3,4,5],4) .== [1.5,2.0,2.5,3.5,4.0])
    @test all(ExponentialNoiseSDE.movmean([1,2,3,4],3) .== [1.5,2.0,3.0,3.5])
    @test all(ExponentialNoiseSDE.movmean([10,2,3,2,4,8,21],4) .== [6.0,5.0,4.25,2.75,4.25,8.75,11.0])
    @test all(ExponentialNoiseSDE.movmean([9,2,4,3,2,7,21],5) .== [5.0,4.5,4.0,3.6,7.4,8.25,10.0])
end

## Small data test
small_data_test()

## Validate calculations (against MATLAB calculations)

X_data = readdlm("./data/exampleData.txt")[:]

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

