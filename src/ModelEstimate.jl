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
    fEstimate, gEstimate, fInitial, gInitial = fg_interate(
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

    # This can be done better with multiple dispatch
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

function theta_search(X,dt,nuMax,thetaMax,betaConvergenceValue)
    # Autocorrelation (single data)
    dA = autocorr_increment(X,nuMax)

    # First search
    thetaStarInitial, rNuMatrix, _ = theta_basis_function_fit(
        dA,
        dt,
        nuMax,
        thetaMax,
        betaConv
    )

    # New maximum tau
    newNuMax = ceil(sqrt(thetaStarInitial)/dt)

    if newNuMax > nuMax
        newNuMax = nuMax
        dA = autocorr_increment(X,newNuMax)
    end

    # Second search
    newThetaMax = newNuMax*dt
    thetaStarNew, _, lambdaStar = theta_basis_function_fit(
        dA[1:newNuMax],
        dt,
        newNuMax,
        newThetaMax,
        betaConv
    )

    # R matrix
    rMatrix = rNuMatrix(thetaStarNew) # Returns r[tau,theta] matrix

    return ThetaProperties(
        nuMax,
        rMatrix,
        dA,
        lambdaStar
    )
end

function theta_basis_function_fit(dA,dt,nuMax,thetaMax,betaConvergenceValue)
    # Objective function from matrix
    rNuMatrix = form_r_matrix(dt,nuMax)
    lambdaFunc = theta_c -> rNuMatrix(theta_c) \ dA
    funcValue = theta_c -> sum((dA .- (rNuMatrix(theta_c) * lambdaFunc(theta_c))).^2)

    # Line search using golden section search and parabolic interpolation
    # Original MATLAB: thetaStar = fminbnd(funcValue,0,thetaMax);
    objective_function(x,g) = funcValue(x[1])
    opt = Opt(:LD_MMA, 1)
    opt.lower_bounds = 0.0
    opt.upper_bounds = thetaMax
    opt.min_objective = objective_function
    funcMin, xmin, retval = optimize(opt, [thetaMax])

    # Best fit lambda vector
    lambdaStar = lambdaFunc(thetaStar);

    return thetaStar,rNuMatrix,lambdaStar
end

function form_r_matrix(dt,nuMax)
    # Basis functions
    r1(tau,theta) = tau - theta*(1 - exp(-tau/theta))
    r2(tau,theta) = tau.^2/2 - theta*r1(tau,theta)
    r3(tau,theta) = tau.^3/6 - theta*r2(tau,theta)
    rArray(tau,theta) = [r1(tau,theta),r2(tau,theta),r3(tau,theta)] # ?
    
    # Functions r matrix (reducing dependencies, new method)
    nu = 1:nuMax
    tau_nu = nu*dt
    rNuMatrix = theta -> rArray(tau_nu,theta)
    return rNuMatrix
end

function fg_interate(lambda1_1,lambda2_1,theta,Xcentre,betaConvergenceValue)
    count_max = 2
    
    # Starting values
    fInitial = lambda1_1
    gInitial = sqrt(abs(lambda2_1))
    fNew = fInitial
    gNew = gInitial
    
    # Iterate until converged
    count = 0
    totalError = Inf
    while (totalError > betaConvergenceValue) && (count < count_max)
        count = count + 1
        
        fOld = fNew
        gOld = gNew
        
        # Find updated values
        [fNew,gNew] = fixed_point_iterate(
            lambda1_1,
            lambda2_1,
            fOld,
            gOld,
            theta,
            Xcentre
        )
        
        totalError = sum((fNew .- fOld).^2) + sum((gNew .- gOld).^2)
    end
   return fNew, gNew, fInitial, gInitial 
end

function lambda_search_linear(moment1,moment2,rMatrix)
    lambda1 = rMatrix \ moment1
    lambda2 = rMatrix \ moment2
    return lambda1, lambda2
end

function fixed_point_iterate(lambda1_1,lambda2_1,f,g,theta,xEvalPoints)
    fGrad = fdiffNU(xEvalPoints,f)
    gGrad = fdiffNU(xEvalPoints,g)

    # Updated functions
    fNew = lambda1_1 .- 0.5*g.*gGrad .- 
                0.5*theta*(fGrad.*g.*gGrad .- f.*gGrad.^2)
    gNew = sqrt(abs(lambda2_1 .- 
                theta*(fGrad.*g.^2 .- f.*g.*gGrad)))
    return fNew, gNew
end

#NOTE: Non-uniform finite-differences (can probably specialize & dispatch)
function fdiffNU(x,F)
    n = length(x)
    m = n-2

    h0 = x[2:n-1] - x[1:n-2] # Backwards steps
    h1 = x[3:n] - x[2:n-1] # Forward steps
    den = h0 .* h1 .* (h0 .+ h1) # Denominator
    h0_2 = h0.^2 # Backwards steps squared
    h1_2 = h1.^2 # Forward steps squared

    # Boundaries
    dFdx[n] = ((h0_2[m] + 2*h0[m]*h1[m]) * F(n) - 
                (h0[m] + h1[m])^2 * F[n-1] + h1_2[m]*F[n-2]) / den[m]
    dFdx[1] = (-(h1_2[1] + 2*h0[1]*h1[1] )* F[1] + 
                (h0[1] + h1[1])^2 * F[2] - h0_2[1]*F[3]) / den[1]

    # Main Body
    dFdx[2:n-1] = (-h1_2.*F[1:n-2] + (h1_2 .- h0_2).*F[2:n-1] + h0_2.*F[3:n]) ./ den;
end