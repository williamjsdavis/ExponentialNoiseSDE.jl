## Test utility functions

# Test utility functions and Statistics
function test_utility_statistics()
    @testset "Statistics and utility functions" begin
        @test all(ExponentialNoiseSDE.movmean([1,2,3,4,5],3) .== [1.5,2.0,3.0,4.0,4.5])
        @test all(ExponentialNoiseSDE.movmean([1,2,3,4,5],2) .== [1.0,1.5,2.5,3.5,4.5])
        @test all(ExponentialNoiseSDE.movmean([1,2,3,4,5],4) .== [1.5,2.0,2.5,3.5,4.0])
        @test all(ExponentialNoiseSDE.movmean([1,2,3,4],3) .== [1.5,2.0,3.0,3.5])
        @test all(ExponentialNoiseSDE.movmean([10,2,3,2,4,8,21],4) .== [6.0,5.0,4.25,2.75,4.25,8.75,11.0])
        @test all(ExponentialNoiseSDE.movmean([9,2,4,3,2,7,21],5) .== [5.0,4.5,4.0,3.6,7.4,8.25,10.0])
    end
end

# Test with a small amount of data
function small_data_test()
    ## Observations
    X_data = DelimitedFiles.readdlm("./data/exampleData.txt")[1:1000]
    N_data = length(X_data)
    obs = Observation(
        X_data,
        0.0050
    )
    @testset "Observation attributes" begin
        @test obs.dt == 0.0050
        @test obs.N == N_data
        @test size(obs.X) == (N_data,)
    end

    ## Conditional moment settings
    nTimeShiftSamplePoints = 6
    nEvalPoints = 5
    xEvalLims = (-0.1,0.1)
    kernel = "Epanechnikov"
    bandwidth = 0.01

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
        @test momentSettings.bandwidth == 0.01
    end

    # Conditional moments
    conditionalMoments = build_moments(obs, momentSettings)

    @testset "ConditionalMoments attributes" begin
        @test size(conditionalMoments.moment1) == momentSize
        @test size(conditionalMoments.moment2) == momentSize
    end
    @testset "ConditionalMoments values" begin
        @test !any(isnan.(conditionalMoments.moment1))
        @test !any(isnan.(conditionalMoments.moment2))
    end

    ## Model estimates

    # Settings
    modelSettings1 = ModelEstimateSettings(displayOutputFlag=false)
    modelSettings2 = ModelEstimateSettings(thetaConvergenceValue=1E-5,displayOutputFlag=true)
    @testset "ModelEstimateSettings" begin
        @test modelSettings1.thetaConvergenceValue == 1E-2
        @test modelSettings2.thetaConvergenceValue == 1E-5
        @test modelSettings1.displayOutputFlag == false
        @test modelSettings2.displayOutputFlag == true
    end

    # Estimating
    modelEstimate = estimate_model(conditionalMoments, modelSettings1)

    @testset "ModelEstimate attributes" begin
        @test size(modelEstimate.driftEstimate) == (nEvalPoints,)
        @test size(modelEstimate.noiseEstimate) == (nEvalPoints,)
        @test size(modelEstimate.driftInitial) == (nEvalPoints,)
        @test size(modelEstimate.noiseInitial) == (nEvalPoints,)
    end
    @testset "ModelEstimate values" begin
        @test modelEstimate.correlationEstimate > -1.0
        @test !any(isnan.(modelEstimate.driftInitial))
        @test !any(isnan.(modelEstimate.noiseInitial))
        @test !any(isnan.(modelEstimate.driftEstimate))
        @test !any(isnan.(modelEstimate.noiseEstimate))
    end

    ## Bootstrapping uncertainties

    # Settings
    bootstrapSettings1 = BootstrapSettings(blockLength=100, nSamples=5, displayOutputFlag=false)
    bootstrapSettings2 = BootstrapSettings()
    @testset "BootstrapSettings" begin
        @test bootstrapSettings1.blockLength == 100
        @test bootstrapSettings1.nSamples == 5
        @test bootstrapSettings1.displayOutputFlag == false
        @test bootstrapSettings2.blockLength == 500
        @test bootstrapSettings2.nSamples == 20
        @test bootstrapSettings2.displayOutputFlag == true
    end

    # Estimating
    bootstrapStatistics = estimate_bootstrap_statistics(modelEstimate,bootstrapSettings1)
    @testset "BootstrapSettings" begin
        @test size(bootstrapStatistics.correlationEstimate[:mean]) == ()
        @test size(bootstrapStatistics.correlationEstimate[:percentiles95]) == (2,)
        @test size(bootstrapStatistics.modelError[:mean]) == ()
        @test size(bootstrapStatistics.modelError[:percentiles95]) == (2,)
        @test size(bootstrapStatistics.driftEstimate[:mean]) == (nEvalPoints,)
        @test size(bootstrapStatistics.driftEstimate[:percentiles95]) == (nEvalPoints,)
        @test size(bootstrapStatistics.driftEstimate[:percentiles95][1]) == (2,)
    end
end

# Get estimated model for a large amount of data (using default settings)
function get_default_estimation_settings()
    ## Conditional moment settings
    nTimeShiftSamplePoints = 15
    nEvalPoints = 20
    xEvalLims = (-1.0,1.0)
    kernel = "Epanechnikov"
    bandwidth = 0.1
    
    # Variables/attributes
    timeShiftSamplePoints = collect(1:nTimeShiftSamplePoints)
    
    momentSettings = ConditionalMomentSettings(
        timeShiftSamplePoints,
        nEvalPoints,
        xEvalLims,
        kernel,
        bandwidth,
    )
    
    ## Model estimates
    # Settings
    modelSettings = ModelEstimateSettings()

    return momentSettings, modelSettings
end

# Get estimated model for a large amount of data (using default settings)
function estimate_large_model()
    # Get settings
    momentSettings, modelSettings = get_default_estimation_settings()

    ## Data
    X_data = readdlm("./data/exampleData.txt")[:]

    obs = Observation(
        X_data,
        0.0050
    )
    
    # Conditional moments
    conditionalMoments = build_moments(obs, momentSettings)
    
    ## Model estimates
    modelEstimate = estimate_model(conditionalMoments, modelSettings)

    return modelEstimate
end

# Get estimated model for a large amount of data (unsmoothed)
function estimate_large_model_unsmoothed()
    # Get settings
    momentSettings, _ = get_default_estimation_settings()
    fitSettings = ModelEstimateSettings(
        functionIterateMethod = FiniteDiff()
    )

    ## Data
    X_data = readdlm("./data/exampleData.txt")[:]

    obs = Observation(
        X_data,
        0.0050
    )

    # Do full estimate
    modelEstimate = estimate_model(obs, momentSettings, fitSettings)
    return modelEstimate
end
