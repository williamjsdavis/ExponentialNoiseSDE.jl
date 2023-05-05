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
#    nCounts::Array{Int64,2}
    moment1::Array{Float64,2}
    moment2::Array{Float64,2}
    xEvalPoints::Array{Float64,1}
    obervation::Observation
    momentSettings
end

# MATLAB equivalent: buildMoments.m
# Only Epan kernel so far
function build_moments(observation::Observation, momentSettings::ConditionalMomentSettings)
    xEvalPoints = calculate_xEvalPoints(momentSettings.evalLims, momentSettings.nEvalPoints)

    kernel = Epaneknikov()
    moment1, moment2 = kbr_moments(
        observation.X, 
        momentSettings.timeShiftSamplePoints, 
        xEvalPoints, 
        momentSettings.bandwidth, 
        kernel
    )

    xEvalPointsCollect = xEvalPoints |> collect
    return ConditionalMoments(
        moment1,
        moment2,
        xEvalPointsCollect,
        observation,
        momentSettings
    )
end
