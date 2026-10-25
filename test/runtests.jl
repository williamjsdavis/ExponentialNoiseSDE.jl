using ExponentialNoiseSDE
using DelimitedFiles
using StableRNGs
using Aqua
using Logging
using Statistics
using Test

# Test data lives next to this file, so tests run from any working directory
datafile(name) = joinpath(@__DIR__, "data", name)

# Run f showing only errors, for code whose expected warnings are not under test
quietly(f) = with_logger(f, ConsoleLogger(stderr, Logging.Error))

# Fixed seed for bootstrap tests (StableRNG streams are stable across Julia versions)
const BOOTSTRAP_SEED = 1

include("test_utils.jl")
include("validate_calculations.jl")
include("simulate.jl")
include("helpers.jl")
include("scientific.jl")

@testset "ExponentialNoiseSDE.jl" begin
    ## Utilities and statistics
    test_utility_statistics()
    test_time_shift_validation()
    test_kernels()
    test_points_without_data()
    test_function_iteration()
    test_bootstrap_statistics()
    test_input_validation()
    test_settings_strategies()

    ## Small dataset
    small_data_test()

    ## Validate full data calculations
    # (against MATLAB smoothed finite diff calculations)
    # One pinned iteration, as in MATLAB, so the iteration has not converged
    modelEstimate = @test_logs (:warn, r"did not converge") estimate_large_model()

    validate_basic(modelEstimate)
    validate_theta(modelEstimate)
    validate_lambda(modelEstimate)
    validate_functions(modelEstimate)

    # Against unsmoothed finite diff calculations
    modelEstimateUnsmoothed = @test_logs (:warn, r"did not converge") estimate_large_model_unsmoothed()
    validate_unsmoothed_functions(modelEstimateUnsmoothed,modelEstimate)

    ## Ground truth from simulated data
    @testset "Ground truth" begin
        scientific_tests()
    end

    ## Package quality
    @testset "Aqua" begin
        Aqua.test_all(ExponentialNoiseSDE)
    end
end
