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
        @test ExponentialNoiseSDE.true_runs([true,true,false,true]) == [1:2, 4:4]
        @test ExponentialNoiseSDE.true_runs([false,true,true,true,false,false]) == [2:4]
        @test ExponentialNoiseSDE.true_runs([false,false]) == UnitRange{Int}[]
    end
    @testset "Autocorrelation" begin
        # Linear (not circular) autocorrelation, by direct summation
        X = randn(StableRNG(2), 50)
        Xdemean = X .- sum(X)/length(X)
        lags = 10
        acfDirect = [sum(Xdemean[1:end-k] .* Xdemean[1+k:end]) for k in 0:lags] ./ sum(abs2, Xdemean)
        @test ExponentialNoiseSDE.myAutocorr(X, lags) ≈ acfDirect rtol=1E-12
    end
end

# A second kernel, to check the kernel setting is the one used
struct BoxKernel <: ExponentialNoiseSDE.Kernel end
(::BoxKernel)(x) = abs(x) < 1 ? 0.5 : 0.0

# Test kernel functions and kernel selection
function test_kernels()
    @testset "Kernels" begin
        # Epanechnikov: normalised, unit variance, support [-sqrt(5), sqrt(5)]
        k = Epanechnikov()
        x = range(-3, 3, length=60001)
        @test sum(k.(x))*step(x) ≈ 1 rtol=1E-6
        @test sum(x.^2 .* k.(x))*step(x) ≈ 1 rtol=1E-6
        @test k(sqrt(5) + 1E-9) == 0.0

        # The kernel in the settings is used
        X = DelimitedFiles.readdlm(datafile("exampleData.txt"))[1:5000]
        obs = Observation(X, 0.0050)
        settings(kernel) = ConditionalMomentSettings(collect(1:6), 5, (-0.5,0.5), kernel, 0.1)
        epanMoments = build_moments(obs, settings(Epanechnikov()))
        boxMoments = build_moments(obs, settings(BoxKernel()))
        @test epanMoments.moment1 != boxMoments.moment1
        # With the box kernel, moment 1 at x = 0 for a shift of 1 is the plain
        # mean of the increments that start within one bandwidth of 0
        near = findall(i -> abs(X[i]) < 0.1, 1:length(X)-1)
        @test boxMoments.moment1[1,3] ≈ sum(X[near .+ 1] .- X[near])/length(near)

        @test_throws ArgumentError settings("Epanechnikov")
    end
end

# Test validation of the time shift lags
function test_time_shift_validation()
    @testset "Time shift validation" begin
        settings(lags) = ConditionalMomentSettings(lags, 5, (-1.0,1.0), Epanechnikov(), 0.1)
        @test settings(collect(1:4)).timeShiftSamplePoints == 1:4
        @test settings(2:2:30).timeShiftSamplePoints == 2:2:30
        @test_throws ArgumentError settings(collect(15:-1:1)) # Not increasing
        @test_throws ArgumentError settings([1,1,2,3,4]) # Repeated
        @test_throws ArgumentError settings([0,1,2,3]) # Not positive
        @test_throws ArgumentError settings([1,2,3]) # Too few to identify theta
    end
end

# Test with a small amount of data
function small_data_test()
    ## Observations
    X_data = DelimitedFiles.readdlm(datafile("exampleData.txt"))[1:1000]
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
    kernel = Epanechnikov()
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
    bootstrapSettings1 = BootstrapSettings(blockLength=100, nSamples=7, displayOutputFlag=false)
    bootstrapSettings2 = BootstrapSettings()
    @testset "BootstrapSettings" begin
        @test bootstrapSettings1.blockLength == 100
        @test bootstrapSettings1.nSamples == 7
        @test bootstrapSettings1.displayOutputFlag == false
        @test bootstrapSettings2.blockLength == 500
        @test bootstrapSettings2.nSamples == 20
        @test bootstrapSettings2.displayOutputFlag == true
    end

    # Estimating. Short resamples of this small dataset often give negative
    # squared noise estimates, which are warned about and set to 0.
    bootstrapStatistics = quietly() do
        estimate_bootstrap_statistics(modelEstimate,bootstrapSettings1; rng=StableRNG(BOOTSTRAP_SEED))
    end
    @testset "Bootstrap reproducibility" begin
        repeatStatistics = quietly() do
            estimate_bootstrap_statistics(modelEstimate,bootstrapSettings1; rng=StableRNG(BOOTSTRAP_SEED))
        end
        @test repeatStatistics.correlationEstimate == bootstrapStatistics.correlationEstimate
        @test repeatStatistics.driftEstimate == bootstrapStatistics.driftEstimate
        @test repeatStatistics.noiseEstimate == bootstrapStatistics.noiseEstimate
    end
    @testset "BootstrapStatistics" begin
        @test size(bootstrapStatistics.correlationEstimate[:mean]) == ()
        @test size(bootstrapStatistics.correlationEstimate[:percentiles95]) == (2,)
        @test size(bootstrapStatistics.modelError[:mean]) == ()
        @test size(bootstrapStatistics.modelError[:percentiles95]) == (2,)
        @test size(bootstrapStatistics.driftEstimate[:mean]) == (nEvalPoints,)
        @test size(bootstrapStatistics.driftEstimate[:percentiles95]) == (nEvalPoints,)
        @test size(bootstrapStatistics.driftEstimate[:percentiles95][1]) == (2,)
    end
    @testset "Bootstrap over many seeds" begin
        # Short resampled series can give a correlation time of exactly 0, which
        # used to crash the estimation. The wider bandwidth keeps data under
        # every evaluation point in every resample.
        wideMomentSettings = ConditionalMomentSettings(
            timeShiftSamplePoints, nEvalPoints, xEvalLims, kernel, 0.05)
        wideModelEstimate = estimate_model(build_moments(obs, wideMomentSettings), modelSettings1)
        nCompleted = quietly() do
            count(1:100) do seed
                estimate_bootstrap_statistics(wideModelEstimate,bootstrapSettings1; rng=StableRNG(seed))
                true
            end
        end
        @test nCompleted == 100

        # With the narrow bandwidth, some resamples have no data under an
        # evaluation point; the bootstrap still completes
        nCompletedSparse = @test_logs (:warn, r"No data within the kernel support") match_mode=:any begin
            count(1:100) do seed
                estimate_bootstrap_statistics(modelEstimate,bootstrapSettings1; rng=StableRNG(seed))
                true
            end
        end
        @test nCompletedSparse == 100
    end
end

# Test evaluation points with no data within the kernel support
function test_points_without_data()
    @testset "Evaluation points without data" begin
        X = readdlm(datafile("exampleData.txt"))[:]
        obs = Observation(X, 0.0050)
        settings(nEvalPoints, evalLims) = ConditionalMomentSettings(
            collect(1:15), nEvalPoints, evalLims, Epanechnikov(), 0.1)

        # Limits wider than the data (-1.55 to 3.41)
        conditionalMoments = @test_logs (:warn, r"No data within the kernel support") build_moments(
            obs, settings(20, (-3.0, 3.0)))
        x = conditionalMoments.xEvalPoints
        noData = [!any(abs.(X[1:end-1] .- xj) .< 0.1*sqrt(5)) for xj in x]
        @test any(noData)

        # NaN only where there is no data
        modelEstimate = estimate_model(conditionalMoments, ModelEstimateSettings())
        @test isnan.(modelEstimate.driftEstimate) == noData
        @test isnan.(modelEstimate.noiseEstimate) == noData
        @test isfinite(modelEstimate.modelError.meanAbsoluteError)

        # Elsewhere, the same as estimating on only the points with data
        firstWithData = findfirst(!, noData)
        @test all(.!noData[firstWithData:end])
        subGrid = estimate_model(obs,
            settings(length(x) - firstWithData + 1, (x[firstWithData], 3.0)),
            ModelEstimateSettings())
        @test modelEstimate.driftEstimate[firstWithData:end] ≈ subGrid.driftEstimate rtol=1E-10
        @test modelEstimate.noiseEstimate[firstWithData:end] ≈ subGrid.noiseEstimate rtol=1E-10

        # The bootstrap completes, with NaN percentiles only where there is no data
        bootstrapStatistics = @test_logs (:warn, r"No data") match_mode=:any estimate_bootstrap_statistics(
            modelEstimate, BootstrapSettings(nSamples=5, displayOutputFlag=false); rng=StableRNG(BOOTSTRAP_SEED))
        driftPercentiles = bootstrapStatistics.driftEstimate[:percentiles95]
        @test all(all.(isnan, driftPercentiles[noData]))
        @test all(all.(isfinite, driftPercentiles[.!noData]))
    end
end

# Get estimated model for a large amount of data (using default settings)
function get_default_estimation_settings()
    ## Conditional moment settings
    nTimeShiftSamplePoints = 15
    nEvalPoints = 20
    xEvalLims = (-1.0,1.0)
    kernel = Epanechnikov()
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
    # Settings: MATLAB's fgIter runs exactly one fixed-point iteration, so
    # comparisons against MATLAB pin the iteration count rather than rely on
    # the stopping rule
    modelSettings = ModelEstimateSettings(
        functionIterateMethod = SmoothedFiniteDiff(functionIterMax = 1)
    )

    return momentSettings, modelSettings
end

# Get estimated model for a large amount of data (MATLAB-equivalent settings)
function estimate_large_model()
    # Get settings
    momentSettings, modelSettings = get_default_estimation_settings()

    ## Data
    X_data = readdlm(datafile("exampleData.txt"))[:]

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
        functionIterateMethod = FiniteDiff(functionIterMax = 1)
    )

    ## Data
    X_data = readdlm(datafile("exampleData.txt"))[:]

    obs = Observation(
        X_data,
        0.0050
    )

    # Do full estimate
    modelEstimate = estimate_model(obs, momentSettings, fitSettings)
    return modelEstimate
end

# Test the drift and noise fixed-point iteration
function test_function_iteration()
    @testset "Drift and noise iteration" begin
        X = readdlm(datafile("exampleData.txt"))[:]
        conditionalMoments = build_moments(Observation(X, 0.0050),
            ConditionalMomentSettings(collect(1:15), 20, (-1.0,1.0), Epanechnikov(), 0.1))
        fit(method) = estimate_model(conditionalMoments, ModelEstimateSettings(functionIterateMethod=method))

        # FiniteDiff is the unsmoothed case of SmoothedFiniteDiff
        @test FiniteDiff <: ExponentialNoiseSDE.FunctionIterateMethod
        @test SmoothedFiniteDiff <: ExponentialNoiseSDE.FunctionIterateMethod
        unsmoothed = fit(FiniteDiff())
        windowOne = fit(SmoothedFiniteDiff(nPointSmooth=1))
        @test unsmoothed.driftEstimate == windowOne.driftEstimate
        @test unsmoothed.noiseEstimate == windowOne.noiseEstimate

        # Default settings converge, in more than one iteration
        for model in (unsmoothed, fit(SmoothedFiniteDiff()))
            @test model.functionDiagnostics.converged
            @test 1 < model.functionDiagnostics.iterations < 20
            @test model.functionDiagnostics.negativeNoiseCount == 0
        end

        # The smoothed iteration reaches small tolerances
        tight = fit(SmoothedFiniteDiff(functionConvergenceValue=1E-10, functionIterMax=50))
        @test tight.functionDiagnostics.converged

        # Not converging is reported
        oneIteration = @test_logs (:warn, r"did not converge") fit(SmoothedFiniteDiff(functionIterMax=1))
        @test !oneIteration.functionDiagnostics.converged
        @test oneIteration.functionDiagnostics.iterations == 1

        # Negative squared noise is set to 0, with a warning
        x = collect(range(-1, 1, length=7))
        lambda2 = [1.0, 1.0, 1.0, -0.5, 1.0, 1.0, 1.0]
        _, g, _, gInitial, diagnostics = @test_logs (:warn, r"Negative squared noise") match_mode=:any begin
            ExponentialNoiseSDE.fg_solve_with_data(zeros(7), lambda2, 0.0, x, FiniteDiff())
        end
        @test gInitial[4] == 0
        @test g[4] == 0
        @test diagnostics.negativeNoiseCount > 0

        # Invalid settings
        @test_throws ArgumentError SmoothedFiniteDiff(nPointSmooth=0)
        @test_throws ArgumentError FiniteDiff(functionIterMax=0)
        @test_throws ArgumentError FiniteDiff(functionConvergenceValue=-1.0)
    end
end

# Test bootstrap summary statistics on hand-made samples
function test_bootstrap_statistics()
    @testset "Bootstrap statistics" begin
        # A model estimate, for the evaluation grid and original estimates
        X = readdlm(datafile("exampleData.txt"))[1:5000]
        modelEstimate = estimate_model(Observation(X, 0.0050),
            ConditionalMomentSettings(collect(1:6), 3, (-0.2,0.2), Epanechnikov(), 0.1),
            ModelEstimateSettings())

        # Skewed samples, with mean 2.2 and median 0
        sampleValues = [0.0, 0.0, 0.0, 1.0, 10.0]
        modelError = ExponentialNoiseSDE.ModelError(0.0, 0.0, 0.0)
        samples = [ExponentialNoiseSDE.BootstrapModelEstimate(v, fill(v, 3), fill(2v, 3), modelError)
            for v in sampleValues]
        statistics = ExponentialNoiseSDE.calculate_sample_statistics(modelEstimate, samples,
            BootstrapSettings(biasCorrection=false, displayOutputFlag=false))

        @test statistics.correlationEstimate[:median] == 0.0
        @test statistics.correlationEstimate[:mean] ≈ 2.2
        @test statistics.driftEstimate[:median] == zeros(3)
        @test statistics.driftEstimate[:mean] ≈ fill(2.2, 3)
        @test statistics.noiseEstimate[:median] == zeros(3)
        @test statistics.noiseEstimate[:mean] ≈ fill(4.4, 3)
    end
end
