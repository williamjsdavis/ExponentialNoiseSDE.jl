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
        @test length(obs) == N_data
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
    modelSettings1 = ModelEstimateSettings()
    modelSettings2 = ModelEstimateSettings(thetaMethod=FixedTheta(0.01), functionIterateMethod=FiniteDiff())
    @testset "ModelEstimateSettings" begin
        @test modelSettings1.thetaMethod == EstimateTheta()
        @test modelSettings1.functionIterateMethod == SmoothedFiniteDiff()
        @test modelSettings2.thetaMethod == FixedTheta(0.01)
        @test modelSettings2.functionIterateMethod == FiniteDiff()
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
    bootstrapSettings1 = BootstrapSettings(blockLength=100, nSamples=7)
    bootstrapSettings2 = BootstrapSettings()
    @testset "BootstrapSettings" begin
        @test bootstrapSettings1.blockLength == 100
        @test bootstrapSettings1.nSamples == 7
        @test bootstrapSettings2.blockLength == 500
        @test bootstrapSettings2.nSamples == 20
        @test bootstrapSettings2.biasCorrection
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
        @test isequal(repeatStatistics.samples.correlationEstimate, bootstrapStatistics.samples.correlationEstimate)
        @test isequal(repeatStatistics.samples.driftEstimate, bootstrapStatistics.samples.driftEstimate)
        @test isequal(repeatStatistics.samples.noiseEstimate, bootstrapStatistics.samples.noiseEstimate)
    end
    @testset "BootstrapStatistics" begin
        @test bootstrapStatistics.correlationEstimate.mean isa Float64
        @test bootstrapStatistics.modelError.lower95 isa Float64
        @test size(bootstrapStatistics.driftEstimate.mean) == (nEvalPoints,)
        @test size(bootstrapStatistics.noiseEstimate.upper95) == (nEvalPoints,)
        # The resampled estimates are kept, one column per resample
        @test size(bootstrapStatistics.samples.correlationEstimate) == (7,)
        @test size(bootstrapStatistics.samples.driftEstimate) == (nEvalPoints, 7)
        @test bootstrapStatistics.correlationEstimate.mean ≈ sum(bootstrapStatistics.samples.correlationEstimate)/7
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
            modelEstimate, BootstrapSettings(nSamples=5); rng=StableRNG(BOOTSTRAP_SEED))
        driftSummary = bootstrapStatistics.driftEstimate
        @test all(isnan, driftSummary.lower95[noData]) && all(isnan, driftSummary.upper95[noData])
        @test all(isfinite, driftSummary.lower95[.!noData]) && all(isfinite, driftSummary.upper95[.!noData])
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
        # A model estimate, for the original estimates
        X = readdlm(datafile("exampleData.txt"))[1:5000]
        modelEstimate = estimate_model(Observation(X, 0.0050),
            ConditionalMomentSettings(collect(1:6), 3, (-0.2,0.2), Epanechnikov(), 0.1),
            ModelEstimateSettings())

        # Skewed samples, with mean 2.2 and median 0
        v = [0.0, 0.0, 0.0, 1.0, 10.0]
        samples = ExponentialNoiseSDE.BootstrapSamples(v, repeat(v', 3), repeat(2v', 3), v)
        summarize(biasCorrection) = ExponentialNoiseSDE.summarize_bootstrap(
            modelEstimate, samples, BootstrapSettings(biasCorrection=biasCorrection))

        statistics = summarize(false)
        @test statistics.correlationEstimate.median == 0.0
        @test statistics.correlationEstimate.mean ≈ 2.2
        @test statistics.driftEstimate.median == zeros(3)
        @test statistics.driftEstimate.mean ≈ fill(2.2, 3)
        @test statistics.noiseEstimate.median == zeros(3)
        @test statistics.noiseEstimate.mean ≈ fill(4.4, 3)
        @test statistics.correlationEstimate.lower95 == quantile(v, 0.025)
        @test statistics.correlationEstimate.upper95 == quantile(v, 0.975)

        # Bias correction shifts only the interval, centring the samples on the
        # original estimate
        θ = modelEstimate.correlationEstimate
        corrected = summarize(true)
        @test corrected.correlationEstimate.mean == statistics.correlationEstimate.mean
        @test corrected.correlationEstimate.lower95 ≈ quantile(v .- (2.2 - θ), 0.025)
        @test corrected.correlationEstimate.upper95 ≈ quantile(v .- (2.2 - θ), 0.975)
        @test corrected.modelError.lower95 == statistics.modelError.lower95
    end
end

# Test the block bootstrap resampling
function test_bootstrap_resampling()
    @testset "Bootstrap resampling" begin
        # Circular blocks of consecutive indices, cut to the series length
        indices = ExponentialNoiseSDE.block_bootstrap_indices(StableRNG(3), 10, 4)
        @test length(indices) == 10
        @test all(in(1:10), indices)
        for blockStart in (1, 5, 9), k in 1:min(3, 10 - blockStart)
            @test indices[blockStart + k] == mod(indices[blockStart] + k - 1, 10) + 1
        end

        # Warnings from the resample fits are reported once per kind, with a count
        X = readdlm(datafile("exampleData.txt"))[1:1000]
        modelEstimate = estimate_model(Observation(X, 0.0050),
            ConditionalMomentSettings(collect(1:6), 5, (-0.1,0.1), Epanechnikov(), 0.05),
            ModelEstimateSettings())
        logs, _ = Test.collect_test_logs() do
            estimate_bootstrap_statistics(modelEstimate, BootstrapSettings(blockLength=100, nSamples=100);
                rng=StableRNG(BOOTSTRAP_SEED))
        end
        warnings = filter(log -> log.level == Logging.Warn, logs)
        @test !isempty(warnings)
        @test all(log -> occursin(r"^\d+ of 100 bootstrap fits: ", log.message), warnings)
        @test length(warnings) == length(unique(log -> split(log.message, ": ", limit=2)[2][1:20], warnings))
    end
end

# Test that invalid inputs are rejected where they enter
function test_input_validation()
    @testset "Input validation" begin
        X = readdlm(datafile("exampleData.txt"))[1:2000]
        @test length(Observation(X, 0.0050)) == 2000
        @test_throws ArgumentError Observation(X, -0.0050)
        @test_throws ArgumentError Observation(X, 0.0)
        @test_throws ArgumentError Observation(X, Inf)
        @test_throws ArgumentError Observation([X; NaN], 0.0050)

        settings(; lags=collect(1:15), n=20, lims=(-1.0,1.0), h=0.1) =
            ConditionalMomentSettings(lags, n, lims, Epanechnikov(), h)
        @test_throws ArgumentError settings(n=2)
        @test_throws ArgumentError settings(lims=(1.0,-1.0))
        @test_throws ArgumentError settings(lims=(-Inf,1.0))
        @test_throws ArgumentError settings(h=0.0)
        @test_throws ArgumentError settings(h=-0.1)
        @test_throws ArgumentError build_moments(Observation(X[1:10], 0.0050), settings())

        @test_throws ArgumentError BootstrapSettings(blockLength=0)
        @test_throws ArgumentError BootstrapSettings(nSamples=1)
        modelEstimate = quietly() do
            estimate_model(Observation(X, 0.0050), settings(), ModelEstimateSettings())
        end
        @test_throws ArgumentError estimate_bootstrap_statistics(modelEstimate,
            BootstrapSettings(blockLength=2000))

        @test_throws ArgumentError FixedTheta(-0.01)
        @test_throws ArgumentError FixedTheta(Inf)
    end
end

# Test the correlation time strategies and progress logging
function test_settings_strategies()
    @testset "Settings and strategies" begin
        X = readdlm(datafile("exampleData.txt"))[:]
        momentSettings = ConditionalMomentSettings(
            timeShiftSamplePoints=1:15, nEvalPoints=20, evalLims=(-1.0,1.0), bandwidth=0.1)
        @test momentSettings.kernel == Epanechnikov()
        conditionalMoments = build_moments(Observation(X, 0.0050), momentSettings)

        # Fixing theta at its estimate gives the same fit
        estimated = estimate_model(conditionalMoments, ModelEstimateSettings())
        fixed = estimate_model(conditionalMoments,
            ModelEstimateSettings(thetaMethod=FixedTheta(estimated.correlationEstimate)))
        @test fixed.correlationEstimate == estimated.correlationEstimate
        @test fixed.driftEstimate ≈ estimated.driftEstimate rtol=1E-12
        @test fixed.noiseEstimate ≈ estimated.noiseEstimate rtol=1E-12
        fixedOther = estimate_model(conditionalMoments, ModelEstimateSettings(thetaMethod=FixedTheta(0.02)))
        @test fixedOther.correlationEstimate == 0.02

        # Bootstrap progress goes to the logging system, not stdout
        printed = mktemp() do path, io
            redirect_stdout(io) do
                @test_logs (:debug, r"Starting") (:debug, "Bootstrap sample 1 of 2") (:debug, "Bootstrap sample 2 of 2") min_level=Logging.Debug match_mode=:any begin
                    estimate_bootstrap_statistics(estimated, BootstrapSettings(nSamples=2); rng=StableRNG(BOOTSTRAP_SEED))
                end
            end
            flush(io)
            read(path, String)
        end
        @test isempty(printed)
    end
end

# Test accessors and printing
function test_api()
    @testset "Accessors and printing" begin
        X = readdlm(datafile("exampleData.txt"))[:]
        observation = Observation(X, 0.0050)
        conditionalMoments = build_moments(observation, ConditionalMomentSettings(
            timeShiftSamplePoints=1:15, nEvalPoints=20, evalLims=(-1.0,1.0), bandwidth=0.1))
        modelEstimate = estimate_model(conditionalMoments, ModelEstimateSettings())
        bootstrapStatistics = estimate_bootstrap_statistics(modelEstimate, BootstrapSettings(nSamples=3);
            rng=StableRNG(BOOTSTRAP_SEED))

        @test correlation_time(modelEstimate) == modelEstimate.correlationEstimate
        @test drift(modelEstimate) == modelEstimate.driftEstimate
        @test noise(modelEstimate) == modelEstimate.noiseEstimate
        @test eval_points(modelEstimate) == conditionalMoments.xEvalPoints
        @test correlation_time(bootstrapStatistics) === bootstrapStatistics.correlationEstimate
        @test drift(bootstrapStatistics).mean == bootstrapStatistics.driftEstimate.mean
        @test noise(bootstrapStatistics).upper95 == bootstrapStatistics.noiseEstimate.upper95

        # Short summaries rather than every array
        for object in (observation, conditionalMoments, modelEstimate, bootstrapStatistics)
            @test length(sprint(show, object)) < 200
            @test length(sprint(show, MIME"text/plain"(), object)) < 2000
        end
        @test occursin("θ = ", sprint(show, modelEstimate))
        @test occursin("converged", sprint(show, MIME"text/plain"(), modelEstimate))
    end
end

# Test the conditional moments against the naive reference implementation
function test_moments_reference()
    @testset "Moments match the naive implementation" begin
        X = readdlm(datafile("exampleData.txt"))[1:10000]
        for lags in (collect(1:15), [1,3,7,20]), lims in ((-1.0,1.0), (-3.0,3.0)),
                bandwidth in (0.01, 0.1), kernel in (Epanechnikov(), BoxKernel())
            x = range(lims[1], lims[2], length=20)
            moment1, moment2 = ExponentialNoiseSDE.kbr_moments(X, lags, x, bandwidth, kernel)
            reference1, reference2 = reference_kbr_moments(X, lags, x, bandwidth, kernel)
            # Same additions in the same order, so identical (NaN where no data)
            @test isequal(moment1, reference1)
            @test isequal(moment2, reference2)
        end
    end
end
