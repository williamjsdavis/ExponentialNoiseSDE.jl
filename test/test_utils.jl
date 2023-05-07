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
    X_data = DelimitedFiles.readdlm("./data/exampleData.txt")[1:100]
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
    modelSettings1 = ModelEstimateSettings()
    modelSettings2 = ModelEstimateSettings(thetaConvergenceValue = 1E-5)
    @testset "ModelEstimateSettings" begin
        @test modelSettings1.thetaConvergenceValue == 1E-2
        @test modelSettings2.thetaConvergenceValue == 1E-5
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
end

# Get estimated model for a large amount of data
function estimate_large_model()
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

    return modelEstimate
end

# Test basic properties
function validate_basic(modelEstimate)
    xEvalPointsMATLAB = DelimitedFiles.readdlm("./data/evalPointsMATLAB.txt")[:]

    xEvalPointsTol = 1E-6
    xEvalPointsError = abs.(modelEstimate.conditionalMoments.xEvalPoints .- xEvalPointsMATLAB)

    @testset "MATLAB validation: basic properties" begin
        @test all(xEvalPointsError .< xEvalPointsTol)
    end
end

# Test theta properties
function validate_theta(modelEstimate)
    # Sample autocorrelation
    dAMATLAB = DelimitedFiles.readdlm("./data/dAMATLAB.txt")[:]

    dATol = 1E-5
    dAError = abs.(modelEstimate.thetaProperties.dA .- dAMATLAB)

    # Lambda star
    # New error (precomputed dA)  1.036447083713199e-8 4.922411280805505e-7 1.0629501136616426e-5
    #NOTE: Possible large error on these coefficients
    lambdaStarMATLAB = DelimitedFiles.readdlm("./data/lambdaStarMATLAB.txt")[:]

    lambdaStarRTol = 0.25
    lambdaStarError = abs.(modelEstimate.thetaProperties.lambdaStar .- lambdaStarMATLAB)
    lambdaStarRError = abs.(lambdaStarError ./ lambdaStarMATLAB)

    # New nMax (should match exactly)
    newNuMaxMATLAB = DelimitedFiles.readdlm("./data/newNuMaxMATLAB.txt", Int64)[1]

    # rMatrix
    rMatrixMATLAB = DelimitedFiles.readdlm("./data/rMatrixMATLAB.txt",',')

    rMatrixTol = 1E-4
    rMatrixError = abs.(modelEstimate.thetaProperties.rMatrix .- rMatrixMATLAB)

    # Correlation estimate
    # MATLAB value      = 0.0102596398898808
    # Previous estimate = 0.01023577455198391
    # New error (precomputed dA) 2.6393070541586017e-10
    correlationEstimteMATLAB = DelimitedFiles.readdlm("./data/correlationEstimateMATLAB.txt")[1]

    correlationEstimteTol = 1E-4
    correlationEstimteError = abs(modelEstimate.correlationEstimate - correlationEstimteMATLAB)

    @testset "MATLAB validation: theta properties" begin
        @test all(dAError .< dATol)
        @test all(lambdaStarRError .< lambdaStarRTol)
        @test newNuMaxMATLAB == modelEstimate.thetaProperties.nuMax
        @test all(rMatrixError .< rMatrixTol)
        @test correlationEstimteError .< correlationEstimteTol
    end
end

# Test lambda properties
function validate_lambda(modelEstimate)
    # lambda1Star/lambda1Est
    lambda1StarMATLAB = DelimitedFiles.readdlm("./data/lambda1StarMATLAB.txt",',')

    lambda1StarRTol = 0.25
    lambda1StarError = abs.(modelEstimate.lambdaProperties.lambda1Est' .- lambda1StarMATLAB)
    lambda1StarRError = abs.(lambda1StarError ./ lambda1StarMATLAB)

    # lambda2Star/lambda2Est
    lambda2StarMATLAB = DelimitedFiles.readdlm("./data/lambda2StarMATLAB.txt",',')

    lambda2StarRTol = 0.12
    lambda2StarError = abs.(modelEstimate.lambdaProperties.lambda2Est' .- lambda2StarMATLAB)
    lambda2StarRError = abs.(lambda2StarError ./ lambda2StarMATLAB)

    @testset "MATLAB validation: lambda properties" begin
        @test all(lambda1StarRError .< lambda1StarRTol)
        @test all(lambda2StarRError .< lambda2StarRTol)
    end
end

# Test estimated functions
function validate_functions(modelEstimate)
    # Drift function (initial)
    driftInitialMATLAB = DelimitedFiles.readdlm("./data/driftInitialMATLAB.txt")[:]

    driftInitialTol = 1E-2
    driftInitialError = abs.(modelEstimate.driftInitial .- driftInitialMATLAB)

    # Drift function (estimate)
    driftEstimateMATLAB = DelimitedFiles.readdlm("./data/driftEstimateMATLAB.txt")[:]

    driftEstimteTol = 1E-2
    driftEstimateError = abs.(modelEstimate.driftEstimate .- driftEstimateMATLAB)

    # Noise function (initial)
    noiseInitialMATLAB = DelimitedFiles.readdlm("./data/noiseInitialMATLAB.txt")[:]

    noiseInitialTol = 1E-3
    noiseInitialError = abs.(modelEstimate.noiseInitial .- noiseInitialMATLAB)

    # Noise function (estimate)
    noiseEstimateMATLAB = DelimitedFiles.readdlm("./data/noiseEstimateMATLAB.txt")[:]

    noiseEstimteTol = 1E-3
    noiseEstimateError = abs.(modelEstimate.noiseEstimate .- noiseEstimateMATLAB)

    @testset "MATLAB validation: function properties" begin
        @test all(driftInitialError .< driftInitialTol)
        @test all(driftEstimateError .< driftEstimteTol)
        @test all(noiseInitialError .< noiseInitialTol)
        @test all(noiseEstimateError .< noiseEstimteTol)
    end
end