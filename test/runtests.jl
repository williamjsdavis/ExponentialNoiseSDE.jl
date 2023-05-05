using ExponentialNoiseSDE
using Test

@testset "First tests" begin
    @test 1==1
end


## Observations
@testset "Observation" begin
    obs = Observation([0.0,0.0,0.0],0.1)
    @test obs.dt == 0.1
end

## Conditional moment settings
timeShiftSamplePoints = collect(1:10)
nEvalPoints = 5
evalLims = (-1.0,1.0)
kernel = "Epanechnikov"
bandwidth = 0.1
conditionalMomentSettings = ConditionalMomentSettings(
    timeShiftSamplePoints,
    nEvalPoints,
    evalLims,
    kernel,
    bandwidth,
)

@testset "ConditionalMomentSettings" begin
    @test conditionalMomentSettings.bandwidth == 0.1
end

#=
# Conditional moments
nSamplePoints = timeShiftSamplePoints |> length
nCounts = zeros(Int64, nSamplePoints, nEvalPoints)
moment1 = zeros(Float64, nSamplePoints, nEvalPoints)
moment2 = zeros(Float64, nSamplePoints, nEvalPoints)
evalPoints = 
obervation
momentOptions
conditionalMomentSettings = ConditionalMomentSettings(
    nCounts,
    moment1
    moment2
    evalPoints
    obervation
    momentOptions
)

@testset "ConditionalMoments" begin
    @test obs.dt == 0.1
end
=#