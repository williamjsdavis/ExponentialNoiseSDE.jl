## Estimation of theta and other properties

# Theta properties
struct ThetaProperties
    nuMax::Int64
    rMatrix::Array{Float64,2}
    dA::Array{Float64,1}
    lambdaStar::Array{Float64,1}
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
            maximumTheta
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

#NOTE: MATLAB secn-PR follows this with a second search over 1:ceil(sqrt(theta)/dt)
# time shifts. That rule depends on the units of time, so it is not used here.
function theta_search(X,dt,nuMax,thetaMax)
    # Autocorrelation (single data)
    dA = autocorr_increment(X,nuMax)

    # Search over all time shifts
    thetaEst, rNuMatrix, lambdaStar = theta_basis_function_fit(
        dA,
        dt,
        nuMax,
        thetaMax
    )

    # R matrix
    rMatrix = rNuMatrix(thetaEst) # Returns r[tau,theta] matrix

    return thetaEst, ThetaProperties(
        nuMax,
        rMatrix,
        dA,
        lambdaStar
    )
end

#NOTE: Make these anon functions into regular functions?
function theta_basis_function_fit(dA,dt,nuMax,thetaMax)
    # Objective function from matrix
    rNuMatrix = form_r_matrix(dt,nuMax)
    lambdaFunc = theta_c -> rNuMatrix(theta_c) \ dA
    funcValue = theta_c -> sum((dA .- (rNuMatrix(theta_c) * lambdaFunc(theta_c))).^2)

    # Bounded line search, MATLAB original: thetaStar = fminbnd(funcValue,0,thetaMax);
    # Searches over u = theta/thetaMax in [0,1], so the tolerance is independent
    # of the units of time
    objective_function(u,g) = funcValue(u[1]*thetaMax)
    opt = Opt(:LN_COBYLA, 1)
    opt.lower_bounds = 0.0
    opt.upper_bounds = 1.0
    opt.xtol_abs = 1E-12
    opt.maxeval = 10_000
    opt.min_objective = objective_function
    funcMin, uMin, retval = optimize(opt, [1.0])
    if retval ∉ (:SUCCESS, :FTOL_REACHED, :XTOL_REACHED)
        @warn "Correlation time search stopped early" retval
    end
    thetaEst = uMin[1]*thetaMax
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