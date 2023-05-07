# Estimating SDE model (drift and noise functions + correlation time)

# MATLAB equivalent: FitOptionsClass.m
struct ModelEstimateSettings
    thetaConvergenceValue::Float64
    functionIterateMethod
    fixThetaFlag::Bool
    fixThetaValue::Float64
    keepObservation::Bool
    displayOutputFlag::Bool
end

# MATLAB equivalent: SPmodelClass.m
struct ModelEstimate
    correlationEstimate::Float64
    driftEstimate::Array{Float64,1}
    noiseEstimate::Array{Float64,1}
    driftInitial::Array{Float64,1}
    noiseInitial::Array{Float64,1}
    modelError
    conditionalMoments::ConditionalMoments
    thetaProperties
    lambdaProperties
    fitSettings::ModelEstimateSettings
end

# MATLAB equivalent: estimateSPmodel.m
function estimate_model(conditionalMoments::ConditionalMoments, fitSettings::ModelEstimateSettings)
    # Estimate correlation time (theta)
    thetaEst, thetaProperties = estimate_theta(conditionalMoments, fitSettings)

    # Estimate drift and noise properties
    lambdaProperties = lambda_search_linear(
        conditionalMoments.moment1,
        conditionalMoments.moment2,
        thetaProperties.rMatrix
    )
    
    # Estimate drift and noise functions
    lambda1_1 = lambdaProperties.lambda1Est[1,:] # Array of lambda^(1)_1
    lambda2_1 = lambdaProperties.lambda2Est[1,:] # Array of lambda^(2)_1
    fEstimate, gEstimate, fInitial, gInitial = fg_solve(
        lambda1_1,
        lambda2_1,
        thetaEst,
        conditionalMoments.xEvalPoints,
        fitSettings.functionIterateMethod
    )

    # Mean fit error
    modelError = mean_fit_error(
        conditionalMoments,
        thetaProperties,
        lambdaProperties,
        fitSettings.displayOutputFlag
    )

    correlationEstimate = thetaEst
    driftEstimate = fEstimate
    noiseEstimate = gEstimate
    driftInitial = fInitial
    noiseInitial = gInitial
    return ModelEstimate(
        correlationEstimate,
        driftEstimate,
        noiseEstimate,
        driftInitial,
        noiseInitial,
        modelError,
        conditionalMoments,
        thetaProperties,
        lambdaProperties,
        fitSettings
    )
end

# Other call args
function estimate_model(observation,momentSettings,fitSettings)
    # Calculate moments
    conditionalMoments = build_moments(observation, momentSettings)
    
    # Estimate model
    modelEstimate = estimate_model(conditionalMoments, fitSettings)
    return modelEstimate
end
