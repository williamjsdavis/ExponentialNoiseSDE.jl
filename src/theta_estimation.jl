## Estimation of theta and other properties

# Theta properties
struct ThetaProperties
    nuMax
    rMatrix # Array{Float64,2} ?
    dA
    lambdaStar
end

# Estimating theta
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
            fitSettings.thetaConvergenceValue,
            fitSettings.displayOutputFlag
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

function theta_search(X,dt,nuMax,thetaMax,thetaConvergenceValue,displayOutputFlag)
    # Autocorrelation (single data)
    dA = autocorr_increment(X,nuMax)

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
    if newNuMax > nuMax
        if displayOutputFlag
            println("Calculating new ACF")
        end
        newNuMax = nuMax
        dA = autocorr_increment(X,newNuMax)
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
    funcMin, xMin, retval = optimize(opt, [thetaMax])
    thetaEst = xMin[1]
    # Best fit lambda vector
    lambdaStar = lambdaFunc(thetaEst)

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