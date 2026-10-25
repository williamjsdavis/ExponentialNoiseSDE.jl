# Conditional moments of observations

# MATLAB equivalent: MomentOptionsClass.m
#TODO: Add optimal bandwidth selection
struct ConditionalMomentSettings
    timeShiftSamplePoints::Array{Int64,1}
    nEvalPoints::Int64
    evalLims::Tuple{Float64,Float64}
    kernel::String # For now
    bandwidth::Float64
    function ConditionalMomentSettings(timeShiftSamplePoints, nEvalPoints, evalLims, kernel, bandwidth)
        check_time_shifts(timeShiftSamplePoints)
        return new(timeShiftSamplePoints, nEvalPoints, evalLims, kernel, bandwidth)
    end
end

# Time shifts are lags in samples: positive and strictly increasing. At least 4
# are needed, since with 3 lags the 3 basis function coefficients fit exactly
# for every correlation time.
function check_time_shifts(timeShiftSamplePoints)
    length(timeShiftSamplePoints) >= 4 || throw(ArgumentError(
        "at least 4 time shifts are needed, got $(length(timeShiftSamplePoints))"))
    all(>(0), diff(timeShiftSamplePoints)) || throw(ArgumentError(
        "time shifts must be strictly increasing, got $(timeShiftSamplePoints)"))
    first(timeShiftSamplePoints) >= 1 || throw(ArgumentError(
        "time shifts must be positive, got $(timeShiftSamplePoints)"))
    return nothing
end

# MATLAB equivalent: MomentClass.m
#Note: Only first 2 moments
struct ConditionalMoments
#    nCounts::Array{Int64,2}
    moment1::Array{Float64,2}
    moment2::Array{Float64,2}
    xEvalPoints::Array{Float64,1}
    timeShifts::Array{Float64,1} # Lags in units of time, one per moment row
    dA::Array{Float64,1} # Autocorrelation increments at the time shifts
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

    # Autocorrelation increments at the same time shifts, for fitting the
    # correlation time
    lags = momentSettings.timeShiftSamplePoints
    timeShifts = lags .* observation.dt
    dA = autocorr_increment(observation.X, maximum(lags))[lags]

    xEvalPointsCollect = xEvalPoints |> collect
    return ConditionalMoments(
        moment1,
        moment2,
        xEvalPointsCollect,
        timeShifts,
        dA,
        observation,
        momentSettings
    )
end
