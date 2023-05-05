# Conditional moments of observations

# MATLAB equivalent: MomentOptionsClass.m
#TODO: Add optimal bandwidth selection
struct ConditionalMomentSettings
    timeShiftSamplePoints::Array{Int64,1}
    nEvalPoints::Int64
    evalLims::Tuple{Float64,Float64}
    kernel::String # For now
    bandwidth::Float64
end

# MATLAB equivalent: MomentClass.m
#Note: Only first 2 moments
struct ConditionalMoments
    nCounts::Array{Int64,2}
    moment1::Array{Float64,2}
    moment2::Array{Float64,2}
    evalPoints::Array{Float64,1}
    obervation::Observation
    momentOptions
end

# MATLAB equivalent: buildMoments.m
function 