using ExponentialNoiseSDE
using DelimitedFiles
using StableRNGs
using Aqua
using Statistics
using Test

# Test data lives next to this file, so tests run from any working directory
datafile(name) = joinpath(@__DIR__, "data", name)

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

    ## Small dataset
    small_data_test()

    ## Validate full data calculations
    # (against MATLAB smoothed finite diff calculations)
    modelEstimate = estimate_large_model()

    validate_basic(modelEstimate)
    validate_theta(modelEstimate)
    validate_lambda(modelEstimate)
    validate_functions(modelEstimate)

    # Against unsmoothed finite diff calculations
    modelEstimateUnsmoothed = estimate_large_model_unsmoothed()
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
