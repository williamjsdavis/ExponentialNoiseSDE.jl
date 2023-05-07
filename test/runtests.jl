using ExponentialNoiseSDE
using DelimitedFiles
using Test
include("test_utils.jl")

## Utilities and statistics
test_utility_statistics()

## Small dataset
small_data_test()

## Validate full data calculations 
# (against MATLAB calculations)
modelEstimate = estimate_large_model()

validate_basic(modelEstimate)
validate_theta(modelEstimate)
validate_lambda(modelEstimate)
validate_functions(modelEstimate)

