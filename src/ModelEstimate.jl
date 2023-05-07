# Estimating SDE model (drift and noise functions + correlation time)

abstract type FunctionIterateMethod end

struct FiniteDiff <: FunctionIterateMethod
    functionConvergenceValue::Float64
    functionIterMax::Int64
    function FiniteDiff(;
        functionConvergenceValue = 0.2,
        functionIterMax = 20
    ) begin new(
        functionConvergenceValue,
        functionIterMax
    )
    end
    end
end

struct SmoothedFiniteDiff 
    nPointSmooth::Int64
    functionConvergenceValue::Float64
    functionIterMax::Int64
    function SmoothedFiniteDiff(;
        nPointSmooth = 3,
        functionConvergenceValue = 0.2,
        functionIterMax = 20
    ) begin new(
        nPointSmooth,
        functionConvergenceValue,
        functionIterMax
    )
    end
    end
end

# MATLAB equivalent: FitOptionsClass.m
struct ModelEstimateSettings
    thetaConvergenceValue::Float64
    functionIterateMethod::Union{FiniteDiff,SmoothedFiniteDiff}
    fixThetaFlag::Bool
    fixThetaValue::Float64
    keepObservation::Bool
    displayOutputFlag::Bool
    # Default settings
    function ModelEstimateSettings(;
        thetaConvergenceValue = 1E-2,
        functionIterateMethod = SmoothedFiniteDiff(),
        fixThetaFlag = false,
        fixThetaValue = 0.0,
        keepObservation = true,
        displayOutputFlag = true
    ) begin
        return new(
            thetaConvergenceValue,
            functionIterateMethod,
            fixThetaFlag,
            fixThetaValue,
            keepObservation,
            displayOutputFlag
        )
    end
    end
end

struct ModelError
    meanAbsoluteErrorMoment1
    meanAbsoluteErrorMoment2
    meanAbsoluteError
end

# MATLAB equivalent: SPmodelClass.m
struct ModelEstimate
    correlationEstimate::Float64
    driftEstimate::Array{Float64,1}
    noiseEstimate::Array{Float64,1}
    driftInitial::Array{Float64,1}
    noiseInitial::Array{Float64,1}
    modelError::ModelError
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

function mean_fit_error(
    conditionalMoments::ConditionalMoments,
    thetaProperties::ThetaProperties,
    lambdaProperties::LambdaProperties,
    displayOutputFlag)

    # Fit differences
    moment1Difference = conditionalMoments.moment1 - 
        thetaProperties.rMatrix * lambdaProperties.lambda1Est
    moment2Difference = conditionalMoments.moment2 - 
        thetaProperties.rMatrix * lambdaProperties.lambda2Est

    # Fit summed mean absolute errors
    meanAbsoluteErrorMoment1 = mean(abs.(moment1Difference))
    meanAbsoluteErrorMoment2 = mean(abs.(moment2Difference))
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

function estimate_theta(conditionalMoments::ConditionalMoments, fitSettings::ModelEstimateSettings)
    
    maximumLag = maximum(conditionalMoments.momentSettings.timeShiftSamplePoints)
    dt = conditionalMoments.obervation.dt

    # This can be done better with multiple dispatch
    if fitSettings.fixThetaFlag
        # Fixed theta method
        thetaEst, thetaProperties = theta_fixed(
            conditionalMoments.obervation.X,
            dt,
            maximumLag,
            fitSettings
        )
    else
        # Estimating theta
        maximumTheta = maximumLag*dt
        thetaEst, thetaProperties = theta_search(
            conditionalMoments.obervation.X,
            dt,
            maximumLag,
            maximumTheta,
            fitSettings.thetaConvergenceValue
        )
    end

    return thetaEst, thetaProperties
end

function theta_fixed(X,dt,nuMax,fitSettings)
    # Get set theta value
    thetaEst = fitSettings.fixThetaValue

    # Autocorrelation (single data)
    dA = autocorr_increment(X,nuMax)

    # R matrix and lambda vector
    rNuMatrix = form_r_matrix(dt,nuMax) # Function: real -> matrix
    rMatrix = rNuMatrix(thetaEst) # Returns r[tau,theta] matrix
    lambdaStar = rMatrix \ dA

    return thetaEst, ThetaProperties(
        nuMax,
        rMatrix,
        dA,
        lambdaStar
    )
end

function theta_search(X,dt,nuMax,thetaMax,thetaConvergenceValue)
    # Autocorrelation (single data)
    @show nuMax
    dA = autocorr_increment(X,nuMax)
    #dA = load_precalculated_dA()

    # First search
    thetaEstInitial, rNuMatrix, _ = theta_basis_function_fit(
        dA,
        dt,
        nuMax,
        thetaMax,
        thetaConvergenceValue
    )
    # New maximum tau
    newNuMax = ceil(Int64, sqrt(thetaEstInitial)/dt)
    @show newNuMax
    if newNuMax > nuMax
        println("Calculating new ACF")
        newNuMax = nuMax
        dA = autocorr_increment(X,newNuMax)
        #dA = load_precalculated_second_dA()
    end

    # Second search
    newThetaMax = newNuMax*dt
    thetaEstNew, _, lambdaStar = theta_basis_function_fit(
        dA[1:newNuMax],
        dt,
        newNuMax,
        newThetaMax,
        thetaConvergenceValue
    )

    # R matrix
    rMatrix = rNuMatrix(thetaEstNew) # Returns r[tau,theta] matrix

    return thetaEstNew, ThetaProperties(
        nuMax,
        rMatrix,
        dA,
        lambdaStar
    )
end

#NOTE: Make these anon functions into regular functions?
function theta_basis_function_fit(dA,dt,nuMax,thetaMax,thetaConvergenceValue)
    # Objective function from matrix
    rNuMatrix = form_r_matrix(dt,nuMax)
    lambdaFunc = theta_c -> rNuMatrix(theta_c) \ dA
    funcValue = theta_c -> sum((dA .- (rNuMatrix(theta_c) * lambdaFunc(theta_c))).^2)

    # Line search using golden section search and parabolic interpolation
    # Original MATLAB: thetaStar = fminbnd(funcValue,0,thetaMax);
    objective_function(x,g) = funcValue(x[1])
    opt = Opt(:LN_COBYLA, 1)
    opt.lower_bounds = 0.0
    opt.upper_bounds = thetaMax
    opt.xtol_abs = 1E-12
    opt.min_objective = objective_function
    @show funcMin, xMin, retval = optimize(opt, [thetaMax])
    thetaEst = xMin[1]
    # Best fit lambda vector
    @show lambdaStar = lambdaFunc(thetaEst)
    
    @show rNuMatrix(thetaEst)
    @show lambdaFunc(thetaEst)
    @show rNuMatrix(0.010259639889881)
    @show lambdaFunc(0.010259639889881)
    @show dA

    return thetaEst, rNuMatrix, lambdaStar
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
    rNuMatrix = theta -> mapreduce(permutedims, vcat, rArray.(tau_nu,theta))
    return rNuMatrix
end

function fg_solve(lambda1_1,lambda2_1,theta,Xcentre,functionIterateMethod::FiniteDiff)
    # Starting values
    fInitial = lambda1_1
    gInitial = sqrt.(abs.(lambda2_1))
    fNew = fInitial
    gNew = gInitial
    
    # Iterate until converged
    count = 0
    totalError = Inf
    while (totalError > functionIterateMethod.functionConvergenceValue) && (count < functionIterateMethod.functionIterMax)
        count = count + 1
        
        fOld = fNew
        gOld = gNew

        # Estimate gradients (finite differences)
        fGrad = fdiffNU(xEvalPoints,fOld)
        gGrad = fdiffNU(xEvalPoints,gOld)
        
        # Find updated values
        fNew, gNew = fixed_point_iterate(
            lambda1_1,
            lambda2_1,
            fOld,
            gOld,
            fGrad,
            gGrad,
            theta
        )
        
        totalError = sum((fNew .- fOld).^2) + sum((gNew .- gOld).^2)
    end
   return fNew, gNew, fInitial, gInitial 
end

moving_average(vs,n) = [sum(@view vs[i:(i+n-1)])/n for i in 1:(length(vs)-(n-1))]

function fg_solve(lambda1_1,lambda2_1,theta,xEvalPoints,functionIterateMethod::SmoothedFiniteDiff)
    # Starting values
    fInitial = lambda1_1
    gInitial = sqrt.(abs.(lambda2_1))
    fNew = movmean(fInitial,functionIterateMethod.nPointSmooth)
    gNew = movmean(gInitial,functionIterateMethod.nPointSmooth)
    
    # Iterate until converged
    count = 0
    totalError = Inf
    while (totalError > functionIterateMethod.functionConvergenceValue) && (count < functionIterateMethod.functionIterMax)
        count = count + 1
        
        fOld = fNew
        gOld = gNew

        # Estimate gradients (finite differences)
        fOld = movmean(fOld,functionIterateMethod.nPointSmooth)
        gOld = movmean(gOld,functionIterateMethod.nPointSmooth)
        fGrad = fdiffNU(xEvalPoints,fOld)
        gGrad = fdiffNU(xEvalPoints,gOld)
        fGrad = movmean(fGrad,functionIterateMethod.nPointSmooth)
        gGrad = movmean(gGrad,functionIterateMethod.nPointSmooth)
        
        # Find updated values
        fNew, gNew = fixed_point_iterate(
            lambda1_1,
            lambda2_1,
            fOld,
            gOld,
            fGrad,
            gGrad,
            theta
        )
        
        totalError = sum((fNew .- fOld).^2) + sum((gNew .- gOld).^2)
    end
   return fNew, gNew, fInitial, gInitial 
end

function lambda_search_linear(moment1,moment2,rMatrix)
    lambda1 = rMatrix \ moment1
    lambda2 = rMatrix \ moment2
    return LambdaProperties(lambda1, lambda2)
end

function fixed_point_iterate(lambda1_1,lambda2_1,f,g,fGrad,gGrad,theta)
    fNew = lambda1_1 .- 0.5*g.*gGrad .- 
                0.5*theta*(fGrad.*g.*gGrad .- f.*gGrad.^2)
    gNew = sqrt.(abs.(lambda2_1 .- 
                theta*(fGrad.*g.^2 .- f.*g.*gGrad)))
    return fNew, gNew
end

#NOTE: Non-uniform finite-differences (can probably specialize & dispatch)
function fdiffNU(x,F)
    n = length(x)
    dFdx = zeros(n)

    h0 = x[2:n-1] - x[1:n-2] # Backwards steps
    h1 = x[3:n] - x[2:n-1] # Forward steps
    den = h0 .* h1 .* (h0 .+ h1) # Denominator
    h0_2 = h0.^2 # Backwards steps squared
    h1_2 = h1.^2 # Forward steps squared

    # Boundaries
    m = n-2
    dFdx[n] = ((h0_2[m] + 2*h0[m]*h1[m]) * F[n] - 
                (h0[m] + h1[m])^2 * F[n-1] + h1_2[m]*F[n-2]) / den[m]
    dFdx[1] = (-(h1_2[1] + 2*h0[1]*h1[1] )* F[1] + 
                (h0[1] + h1[1])^2 * F[2] - h0_2[1]*F[3]) / den[1]

    # Main Body
    dFdx[2:n-1] = (-h1_2.*F[1:n-2] + (h1_2 .- h0_2).*F[2:n-1] + h0_2.*F[3:n]) ./ den
    
    return dFdx
end

function load_precalculated_dA()
    dA = [
        -0.000582057549633207
        -0.00200082119911433
        -0.00391300658940873
        -0.00611098179579842
        -0.00846928666848296
        -0.0109122835560818
        -0.0133937924967669
        -0.015882948452609
        -0.018358366667826
        -0.0208063744974986
        -0.0232204816798446
        -0.0255967317691245
        -0.027936617102445
        -0.0302447274405242
        -0.0325229560133573
    ]
    return dA
end
function load_precalculated_second_dA()
    dA = [
        -0.000582057549633207
        -0.00200082119911433
        -0.00391300658940873
        -0.00611098179579842
        -0.00846928666848297
        -0.0109122835560818
        -0.0133937924967669
        -0.0158829484526090
        -0.0183583666678260
        -0.0208063744974986
        -0.0232204816798446
        -0.0255967317691245
        -0.0279366171024450
        -0.0302447274405242
        -0.0325229560133573
    ]
    return dA
end