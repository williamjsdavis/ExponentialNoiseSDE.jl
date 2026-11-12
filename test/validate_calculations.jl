## Validate calculations against MATLAB calculations

# Test basic properties
function validate_basic(modelEstimate)
    xEvalPointsMATLAB = DelimitedFiles.readdlm(datafile("evalPointsMATLAB.txt"))[:]

    xEvalPointsTol = 1E-6
    xEvalPointsError = abs.(modelEstimate.conditionalMoments.xEvalPoints .- xEvalPointsMATLAB)

    @testset "MATLAB validation: basic properties" begin
        @test all(xEvalPointsError .< xEvalPointsTol)
    end
end

# MATLAB estimates the autocorrelation increments from the sample
# autocorrelation, and the package from squared increments (#36); the two
# differ by O(lag/N). The other comparisons use MATLAB's autocorrelation
# increments as input (with_matlab_dA), so that they stay exact.
function validate_autocorrelation_increments()
    momentSettings, _ = get_default_estimation_settings()
    observation = Observation(readdlm(datafile("exampleData.txt"))[:], 0.0050)
    dA = build_moments(observation, momentSettings).dA
    dAMATLAB = DelimitedFiles.readdlm(datafile("dAMATLAB.txt"))[:]

    # Agreement achieved: 1.6e-3 relative at lag 1, the largest
    @testset "MATLAB validation: autocorrelation increments" begin
        @test all(abs.(dA ./ dAMATLAB .- 1) .< 5E-3)
    end
end

with_matlab_dA(cm) = ExponentialNoiseSDE.ConditionalMoments(cm.moment1, cm.moment2,
    cm.xEvalPoints, cm.timeShifts, DelimitedFiles.readdlm(datafile("dAMATLAB.txt"))[:],
    cm.observation, cm.momentSettings)

# Test theta properties
function validate_theta(modelEstimate)
    # Tolerances below are 20-600x the agreement achieved (lambdaStar 5.9e-6
    # rel., rMatrix 2.5e-10, theta 2.6e-10)

    # Lambda star
    lambdaStarMATLAB = DelimitedFiles.readdlm(datafile("lambdaStarMATLAB.txt"))[:]

    lambdaStarRTol = 1E-4
    lambdaStarError = abs.(modelEstimate.thetaProperties.lambdaStar .- lambdaStarMATLAB)
    lambdaStarRError = abs.(lambdaStarError ./ lambdaStarMATLAB)

    # New nMax (should match exactly)
    newNuMaxMATLAB = DelimitedFiles.readdlm(datafile("newNuMaxMATLAB.txt"), Int64)[1]

    # rMatrix
    rMatrixMATLAB = DelimitedFiles.readdlm(datafile("rMatrixMATLAB.txt"),',')

    rMatrixTol = 1E-8
    rMatrixError = abs.(modelEstimate.thetaProperties.rMatrix .- rMatrixMATLAB)

    # Correlation estimate
    correlationEstimteMATLAB = DelimitedFiles.readdlm(datafile("correlationEstimateMATLAB.txt"))[1]

    correlationEstimteTol = 1E-8
    correlationEstimteError = abs(modelEstimate.correlationEstimate - correlationEstimteMATLAB)

    @testset "MATLAB validation: theta properties" begin
        @test all(lambdaStarRError .< lambdaStarRTol)
        @test newNuMaxMATLAB == modelEstimate.thetaProperties.nuMax
        @test all(rMatrixError .< rMatrixTol)
        @test correlationEstimteError .< correlationEstimteTol
    end
end

# Test lambda properties
function validate_lambda(modelEstimate)
    # lambda1Star/lambda1Est
    lambda1StarMATLAB = DelimitedFiles.readdlm(datafile("lambda1StarMATLAB.txt"),',')

    # Agreement achieved: 4.9e-6 rel. (lambda1), 2.5e-6 rel. (lambda2)
    lambda1StarRTol = 1E-4
    lambda1StarError = abs.(modelEstimate.lambdaProperties.lambda1Est' .- lambda1StarMATLAB)
    lambda1StarRError = abs.(lambda1StarError ./ lambda1StarMATLAB)

    # lambda2Star/lambda2Est
    lambda2StarMATLAB = DelimitedFiles.readdlm(datafile("lambda2StarMATLAB.txt"),',')

    lambda2StarRTol = 1E-4
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
    driftInitialMATLAB = DelimitedFiles.readdlm(datafile("driftInitialMATLAB.txt"))[:]

    # Agreement achieved: drift 4.1e-8, noise 1.1e-8
    driftInitialTol = 1E-6
    driftInitialError = abs.(modelEstimate.driftInitial .- driftInitialMATLAB)

    # Drift function (estimate)
    driftEstimateMATLAB = DelimitedFiles.readdlm(datafile("driftEstimateMATLAB.txt"))[:]

    driftEstimteTol = 1E-6
    driftEstimateError = abs.(modelEstimate.driftEstimate .- driftEstimateMATLAB)

    # Noise function (initial)
    noiseInitialMATLAB = DelimitedFiles.readdlm(datafile("noiseInitialMATLAB.txt"))[:]

    noiseInitialTol = 1E-6
    noiseInitialError = abs.(modelEstimate.noiseInitial .- noiseInitialMATLAB)

    # Noise function (estimate)
    noiseEstimateMATLAB = DelimitedFiles.readdlm(datafile("noiseEstimateMATLAB.txt"))[:]

    noiseEstimteTol = 1E-6
    noiseEstimateError = abs.(modelEstimate.noiseEstimate .- noiseEstimateMATLAB)

    @testset "MATLAB validation: function properties" begin
        @test all(driftInitialError .< driftInitialTol)
        @test all(driftEstimateError .< driftEstimteTol)
        @test all(noiseInitialError .< noiseInitialTol)
        @test all(noiseEstimateError .< noiseEstimteTol)
    end
end

# Test estimated functions
function validate_unsmoothed_functions(modelEstimateUnsmoothed,modelEstimate)
    # Drift function (estimate)
    driftEstimateMATLAB = DelimitedFiles.readdlm(datafile("driftEstimateNoSmoothMATLAB.txt"))[:]

    # Agreement achieved: drift 4.7e-8, noise 1.2e-8
    driftEstimteTol = 1E-6
    driftEstimateError = abs.(modelEstimateUnsmoothed.driftEstimate .- driftEstimateMATLAB)

    # Noise function (estimate)
    noiseEstimateMATLAB = DelimitedFiles.readdlm(datafile("noiseEstimateNoSmoothMATLAB.txt"))[:]

    noiseEstimteTol = 1E-6
    noiseEstimateError = abs.(modelEstimateUnsmoothed.noiseEstimate .- noiseEstimateMATLAB)

    @testset "MATLAB validation: unsmoothed function properties" begin
        @test all(modelEstimateUnsmoothed.driftInitial .== modelEstimate.driftInitial)
        @test all(modelEstimateUnsmoothed.noiseInitial .== modelEstimate.noiseInitial)
        @test !all(modelEstimateUnsmoothed.driftEstimate .== modelEstimate.driftEstimate)
        @test !all(modelEstimateUnsmoothed.noiseEstimate .== modelEstimate.noiseEstimate)
        @test all(driftEstimateError .< driftEstimteTol)
        @test all(noiseEstimateError .< noiseEstimteTol)
    end
end