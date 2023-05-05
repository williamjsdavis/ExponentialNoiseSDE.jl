# Estimating SDE model (drift and noise functions + correlation time)

# MATLAB equivalent: FitOptionsClass.m
struct ModelEstimateSettings
    thetaConvergenceValue::Float64
    functionConvergenceValue::Float64
    fixThetaFlag::Bool
    fixThetaValue::Float64
    keepObservation::Bool
    displayOutputFlag::Bool
    # Default settings
    function ModelEstimateSettings(;
        thetaConvergenceValue = 1E-2,
        functionConvergenceValue = 0.2,
        fixThetaFlag = false,
        fixThetaValue = 0.0,
        keepObservation = true,
        displayOutputFlag = true
    ) begin
        return new(
            thetaConvergenceValue,
            functionConvergenceValue,
            fixThetaFlag,
            fixThetaValue,
            keepObservation,
            displayOutputFlag
        )
    end
    end
end

# MATLAB equivalent: SPmodelClass.m
struct ModelEstimate
    correlationEstimate::Float64
    driftEstimate::Array{Float64,1}
    noiseEstimate::Array{Float64,1}
    driftInitial::Array{Float64,1}
    noiseInitial::Array{Float64,1}
    meanAbsoluteError::Float64
    conditionalMoments::ConditionalMoments
    thetaProperties
    lambdaProperties
    fitSettings::ModelEstimateSettings
end

struct ThetaProperties
    nuMax
    rMatrix # Array{Float64,2} ?
    dA
    lambdaStar
end
struct LambdaProperties
    lambda1Est # Array{Float64,2} ?
    lambda2Est # Array{Float64,2} ?
end

# MATLAB equivalent: estimateSPmodel.m
function estimate_model(conditionalMoments::ConditionalMoments, fitSettings::ModelEstimateSettings)
    # Estimate correlation time (theta)
    thetaEst, thetaProperties = estimate_theta(conditionalMoments, fitSettings)

    # Estimate drift and noise properties
    lambdaProperties = lambdaSearchLinear(
        conditionalMoments.moment1,
        conditionalMoments.moment2,
        thetaProperties.rMatrix
    )
    
    # Estimate drift and noise functions
    lambda1_1 = lambdaProperties.lambda1Star[1,:] # Array of lambda^(1)_1
    lambda2_1 = lambdaProperties.lambda2Star[1,:] # Array of lambda^(2)_1
    fEstimate, gEstimate, fInitial, gInitial = fgIter(
        lambda1_1,
        lambda2_1,
        thetaEst,
        conditionalMoments.xEvalPoints,
        fitSettings.functionConvergence
    )

    # Mean fit error
    meanAbsoluteError = mean_fit_error(
        conditionalMoments,
        thetaProperties,
        lambdaProperties,
        fitSettings.displayOutputFlag
    )


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
        meanAbsoluteError,
        conditionalMoments,
        thetaProperties,
        lambdaProperties,
        fitOptions
    )
end

struct ModelError
    meanAbsoluteErrorMoment1
    meanAbsoluteErrorMoment2
    meanAbsoluteError
end

function mean_fit_error(
    conditionalMoments::ConditionalMoments,
    thetaProperties::ThetaProperties,
    lambdaProperties::LambdaProperties,
    displayOutputFlag)

    # Fit differences
    moment1Difference = conditionalMoments.moment1 - 
        thetaProperties.rMatrix * lambdaProperties.lambda1Star
    moment2Difference = conditionalMoments.moment2 - 
        thetaProperties.rMatrix * lambdaProperties.lambda2Star

    # Fit summed mean absolute errors
    meanAbsoluteErrorMoment1 = sum(abs(moment1Difference(:))) / length(moment1Difference)
    meanAbsoluteErrorMoment2 = sum(abs(moment2Difference(:))) / length(moment2Difference)
    meanAbsoluteError = 0.5*(meanAbsoluteErrorMoment1 + meanAbsoluteErrorMoment2)

    if displayOutputFlag
        #=
        println('Absolute fit error on moments')
        println(['Mean M^(1) fit error: ',
            sprintf('%0.4e',meanAbsoluteError.moment1)])
        println(['Mean M^(2) fit error: ',
            sprintf('%0.4e',meanAbsoluteError.moment2)])
        println(['Mean fit error: ',
            sprintf('%0.4e',meanAbsoluteError.bothMoments)])
        =#
    end
    
    return ModelError(
        meanAbsoluteErrorMoment1,
        meanAbsoluteErrorMoment2,
        meanAbsoluteError
    )
end

function estimate_theta(conditionalMoments::ConditionalMoments, fitOptions::ModelEstimateSettings)
    
    maximumLag = maximum(conditionalMoments.momentSettings.timeShiftSamplePoints)
    dt = conditionalMoments.obervation.dt


    if fitSettings.fixThetaFlag
        # Fixed theta method
        thetaStar, thetaProperties = theta_fixed(
            conditionalMoments.obervation.X,
            dt,
            maximumLag,
            fitOptions
        )
    else
        # Estimating theta
        maximumTheta = maximumLag*dt
        thetaStar, thetaProperties = theta_search(
            conditionalMoments.obervation.X,
            dt,
            maximumLag,
            maximumTheta,
            fitOptions.thetaConvergence
        )
    end

    return thetaEst, thetaProperties
end

function theta_fixed(X,dt,nuMax,fitSettings)
    # Get set theta value
    thetaStar = fitSettings.fixThetaValue

    # Autocorrelation (single data)
    dA = autocorr_increment(X,nuMax)

    # R matrix and lambda vector
    rNuMatrix = form_r_matrix(dt,nuMax) # Function: real -> matrix
    rMatrix = rNuMatrix(thetaStar) # Returns r[tau,theta] matrix
    lambdaStar = rMatrix \ dA

    return ThetaProperties(
        nuMax,
        rMatrix,
        dA,
        lambdaStar
    )
end

