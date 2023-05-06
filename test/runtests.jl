using ExponentialNoiseSDE
using DelimitedFiles
using Test

@testset "First tests" begin
    @test 1==1
end


## Observations
X_data = DelimitedFiles.readdlm("./data/exampleData.txt")[1:100]
N_data = length(X_data)
obs = Observation(
    X_data,
    0.0050
)
@testset "Observation" begin
    @test obs.dt == 0.0050
end

## Conditional moment settings
nTimeShiftSamplePoints = 6
nEvalPoints = 5
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

@testset "ConditionalMomentSettings" begin
    @test momentSettings.bandwidth == 0.1
end

# Conditional moments
conditionalMoments = build_moments(obs, momentSettings)

@testset "ConditionalMoments" begin
    @test size(conditionalMoments.moment1) == momentSize
end

## Model estimates

# Settings
modelSettings1 = ModelEstimateSettings()
modelSettings2 = ModelEstimateSettings(thetaConvergenceValue = 1E-5)
@testset "ModelEstimateSettings" begin
    @test modelSettings1.thetaConvergenceValue == 1E-2
    @test modelSettings2.thetaConvergenceValue == 1E-5
end

# Estimating
modelEstimate = estimate_model(conditionalMoments, modelSettings1)

@testset "ModelEstimate" begin
    @test size(modelEstimate.driftEstimate) == (nEvalPoints,)
    @test size(modelEstimate.noiseEstimate) == (nEvalPoints,)
    @test size(modelEstimate.driftInitial) == (nEvalPoints,)
    @test size(modelEstimate.noiseInitial) == (nEvalPoints,)
    @test modelEstimate.correlationEstimate > -1.0
end

