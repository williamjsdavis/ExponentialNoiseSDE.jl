## Scientific tests against simulated ground truth
# Tolerances come from test/calibration/calibrate.jl (1.5 × the 99th
# percentile over 100 seeds). Rerun it when the estimator changes.

const θ_TRUE = 0.01
const DT_SIM = 0.005

function scientific_tests()
    # Linear process with constant noise: has an analytic variance and cannot diverge
    a = 2.5
    linear_drift(x) = -a*x
    constant_noise(x) = 1.0

    @testset "Simulator" begin
        X = simulate_process(StableRNG(1), linear_drift, constant_noise, θ_TRUE, DT_SIM, 200_000)
        @test isapprox(var(X), 1/(2a*(1 + a*θ_TRUE)); rtol=0.11)
    end

    # Example model: f = -2.5x, g = 1 + x²/4 for x > 0
    X = simulate_process(StableRNG(1), example_drift, example_noise, θ_TRUE, DT_SIM, 50_000)
    model = fit_model(X, DT_SIM)
    x = eval_points_of(model)
    inner = abs.(x) .<= 0.6

    @testset "Recovery" begin
        @test isapprox(theta_of(model), θ_TRUE; rtol=0.32)
        @test maximum(abs, drift_of(model)[inner] .- example_drift.(x[inner])) < 0.72
        @test maximum(abs, noise_of(model)[inner] .- example_noise.(x[inner])) < 0.15
    end

    @testset "Consistency" begin
        spread(n) = std([theta_of(fit_model(
            simulate_process(StableRNG(s), linear_drift, constant_noise, θ_TRUE, DT_SIM, n), DT_SIM))
            for s in 1:8])
        # Expect ≈ 0.25 for 16× the data; 1 would mean no improvement
        @test spread(160_000)/spread(10_000) < 0.6
    end

    @testset "White noise" begin
        # θ = 0: the estimate should be at or near zero (θ̂/dt ≤ 0.04 over 100
        # seeds; 40% are exactly 0, which used to crash the estimation)
        for seed in 1:10
            Xwhite = simulate_white_noise_process(StableRNG(seed), a, DT_SIM, 50_000)
            @test theta_of(fit_model(Xwhite, DT_SIM)) < 0.1*DT_SIM
        end
    end

    # Estimates must not depend on the units of time or space
    smoothed = SmoothedFiniteDiff()
    unsmoothed = FiniteDiff()
    @testset "Invariance, fixed iterations" begin
        for method in (pinned(smoothed, 3), pinned(unsmoothed, 3))
            reference = fit_model(X, DT_SIM; method)
            @test matches(reference, time_rescaled(X, DT_SIM, method))
            @test matches(reference, space_rescaled(X, DT_SIM, method))
        end
    end
    @testset "Invariance, default settings" begin
        # The default stopping thresholds are absolute, so the number of
        # iterations can change with the units. The smoothed time check passes
        # only because both runs stop after one iteration on this data.
        reference = fit_model(X, DT_SIM; method=smoothed)
        @test matches(reference, time_rescaled(X, DT_SIM, smoothed))
        @test_broken matches(reference, space_rescaled(X, DT_SIM, smoothed))
        reference = fit_model(X, DT_SIM; method=unsmoothed)
        @test_broken matches(reference, time_rescaled(X, DT_SIM, unsmoothed))
        @test_broken matches(reference, space_rescaled(X, DT_SIM, unsmoothed))
    end
end
