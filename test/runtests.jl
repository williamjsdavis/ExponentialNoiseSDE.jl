using ExponentialNoiseSDE
using Test

@testset "First tests" begin
    @test 1==1
end


## Observations
obs = Observation(
    zeros(Float64,100),
    0.1
)
@testset "Observation" begin
    @test obs.dt == 0.1
end

## Conditional moment settings
nTimeShiftSamplePoints = 60
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
estimate_model(conditionalMoments, modelSettings1)

