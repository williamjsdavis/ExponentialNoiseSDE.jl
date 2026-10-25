## Estimation of theta and other properties

# Theta properties
struct ThetaProperties
    nuMax::Int64
    rMatrix::Array{Float64,2}
    dA::Array{Float64,1}
    lambdaStar::Array{Float64,1}
end

# Correlation time (theta) methods
abstract type ThetaMethod end

"""
    EstimateTheta()

Estimate the correlation time θ by a least-squares fit of the basis functions
r_k(τ, θ) to the autocorrelation increments, for θ between 0 and the largest
time shift. The default `thetaMethod`.
"""
struct EstimateTheta <: ThetaMethod end

"""
    FixedTheta(θ)

Use a known correlation time `θ ≥ 0` (in units of time) instead of estimating it.
"""
struct FixedTheta <: ThetaMethod
    theta::Float64
    function FixedTheta(theta)
        (theta >= 0 && isfinite(theta)) || throw(ArgumentError(
            "the fixed correlation time must be non-negative and finite, got $theta"))
        return new(theta)
    end
end

# Estimating theta, from the autocorrelation increments at the moment time shifts
function estimate_theta(conditionalMoments::ConditionalMoments, ::EstimateTheta)
    return theta_search(
        conditionalMoments.dA,
        conditionalMoments.timeShifts,
        maximum(conditionalMoments.momentSettings.timeShiftSamplePoints)
    )
end

function estimate_theta(conditionalMoments::ConditionalMoments, thetaMethod::FixedTheta)
    return theta_fixed(
        conditionalMoments.dA,
        conditionalMoments.timeShifts,
        maximum(conditionalMoments.momentSettings.timeShiftSamplePoints),
        thetaMethod.theta
    )
end

function theta_fixed(dA,timeShifts,nuMax,thetaEst)
    # R matrix and lambda vector
    rMatrix = form_r_matrix(timeShifts,thetaEst) # Returns r[tau,theta] matrix
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
function theta_search(dA,timeShifts,nuMax)
    # Search over all time shifts, for theta up to the largest time shift
    thetaEst, lambdaStar = theta_basis_function_fit(
        dA,
        timeShifts,
        maximum(timeShifts)
    )

    # R matrix
    rMatrix = form_r_matrix(timeShifts,thetaEst) # Returns r[tau,theta] matrix

    return thetaEst, ThetaProperties(
        nuMax,
        rMatrix,
        dA,
        lambdaStar
    )
end

#NOTE: Make these anon functions into regular functions?
function theta_basis_function_fit(dA,timeShifts,thetaMax)
    # Objective function from matrix
    lambdaFunc = theta_c -> form_r_matrix(timeShifts,theta_c) \ dA
    funcValue = theta_c -> sum((dA .- (form_r_matrix(timeShifts,theta_c) * lambdaFunc(theta_c))).^2)

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

    return thetaEst, lambdaStar
end

# Matrix of basis functions r_k(tau,theta), one row per time shift tau
function form_r_matrix(timeShifts,theta)
    rMatrix = Matrix{Float64}(undef, length(timeShifts), 3)
    for (i, tau) in enumerate(timeShifts)
        r1 = tau - theta*(1 - exp(-tau/theta))
        r2 = tau^2/2 - theta*r1
        r3 = tau^3/6 - theta*r2
        rMatrix[i,1] = r1
        rMatrix[i,2] = r2
        rMatrix[i,3] = r3
    end
    return rMatrix
end
