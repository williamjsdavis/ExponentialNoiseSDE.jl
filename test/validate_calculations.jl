## Validate calculations against MATLAB calculations

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

# Test estimated functions
function validate_unsmoothed_functions(modelEstimateUnsmoothed,modelEstimate)
    # Drift function (estimate)
    driftEstimateMATLAB = DelimitedFiles.readdlm("./data/driftEstimateNoSmoothMATLAB.txt")[:]

    driftEstimteTol = 0.15
    driftEstimateError = abs.(modelEstimateUnsmoothed.driftEstimate .- driftEstimateMATLAB)

    # Noise function (estimate)
    noiseEstimateMATLAB = DelimitedFiles.readdlm("./data/noiseEstimateNoSmoothMATLAB.txt")[:]

    noiseEstimteTol = 0.03
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