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

# Test theta properties
function validate_theta(modelEstimate)
    # Sample autocorrelation
    dAMATLAB = DelimitedFiles.readdlm(datafile("dAMATLAB.txt"))[:]

    # Tolerances below are 20-600x the agreement achieved (dA 1.7e-16, lambdaStar
    # 5.9e-6 rel., rMatrix 2.5e-10, theta 2.6e-10)
    dATol = 1E-13
    dAError = abs.(modelEstimate.thetaProperties.dA .- dAMATLAB)

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