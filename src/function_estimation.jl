# Estimation of functions f and g

# Function estimation methods
abstract type FunctionIterateMethod end

# Fixed-point iteration with finite-difference gradients. Stops when the
# relative change in both f and g from the previous iteration is below
# functionConvergenceValue, or after functionIterMax iterations.
struct FiniteDiff <: FunctionIterateMethod
    functionConvergenceValue::Float64
    functionIterMax::Int64
    function FiniteDiff(functionConvergenceValue, functionIterMax)
        check_iteration_settings(functionConvergenceValue, functionIterMax)
        return new(functionConvergenceValue, functionIterMax)
    end
end

# As FiniteDiff, with an nPointSmooth-point moving mean applied to the
# functions and their gradients in each iteration
struct SmoothedFiniteDiff <: FunctionIterateMethod
    nPointSmooth::Int64
    functionConvergenceValue::Float64
    functionIterMax::Int64
    function SmoothedFiniteDiff(nPointSmooth, functionConvergenceValue, functionIterMax)
        nPointSmooth >= 1 || throw(ArgumentError(
            "nPointSmooth must be at least 1, got $nPointSmooth"))
        check_iteration_settings(functionConvergenceValue, functionIterMax)
        return new(nPointSmooth, functionConvergenceValue, functionIterMax)
    end
end

function check_iteration_settings(functionConvergenceValue, functionIterMax)
    functionConvergenceValue >= 0 || throw(ArgumentError(
        "functionConvergenceValue must be non-negative, got $functionConvergenceValue"))
    functionIterMax >= 1 || throw(ArgumentError(
        "functionIterMax must be at least 1, got $functionIterMax"))
    return nothing
end

# Diagnostics of the drift and noise iteration
struct FunctionIterationDiagnostics
    iterations::Int64 # Iterations used (the most over runs of points with data)
    converged::Bool
    negativeNoiseCount::Int64 # Times the squared noise was negative and set to 0
end

# Moving mean window; FiniteDiff is the unsmoothed case
smoothing_window(::FiniteDiff) = 1
smoothing_window(functionIterateMethod::SmoothedFiniteDiff) = functionIterateMethod.nPointSmooth

# Solve separately on each run of consecutive evaluation points with data, so
# points without data (NaN lambdas) do not spread to their neighbours through
# the finite differences and smoothing
function fg_solve_with_data(lambda1_1,lambda2_1,theta,xEvalPoints,functionIterateMethod)
    fInitial = copy(lambda1_1)
    gInitial = sqrt.(max.(lambda2_1, 0))
    fEstimate = fill(NaN, length(xEvalPoints))
    gEstimate = fill(NaN, length(xEvalPoints))
    iterations = 0
    converged = true
    nNegative = count(<(0), lambda2_1)

    hasData = isfinite.(lambda1_1) .& isfinite.(lambda2_1)
    for run in true_runs(hasData)
        length(run) >= 3 || continue # Finite differences need 3 points
        fEstimate[run], gEstimate[run], runIterations, runConverged, runNegative = fg_solve(
            lambda1_1[run],
            lambda2_1[run],
            theta,
            xEvalPoints[run],
            functionIterateMethod
        )
        iterations = max(iterations, runIterations)
        converged &= runConverged
        nNegative += runNegative
    end

    if nNegative > 0
        @warn "Negative squared noise estimate ($nNegative time(s), over points and iterations); set to 0"
    end
    if !converged
        @warn "Drift and noise iteration did not converge in $(functionIterateMethod.functionIterMax) iterations" functionIterateMethod.functionConvergenceValue
    end
    diagnostics = FunctionIterationDiagnostics(iterations, converged, nNegative)
    return fEstimate, gEstimate, fInitial, gInitial, diagnostics
end

# Fixed-point solution for f and g, starting from the lambdas. Returns f, g,
# the number of iterations, whether it converged, and the number of times the
# squared noise was negative (and set to 0).
function fg_solve(lambda1_1,lambda2_1,theta,xEvalPoints,functionIterateMethod::FunctionIterateMethod)
    nPointSmooth = smoothing_window(functionIterateMethod)
    smooth(F) = movmean(F,nPointSmooth)
    tolerance = functionIterateMethod.functionConvergenceValue

    # Starting values
    fNew = smooth(lambda1_1)
    gNew = smooth(sqrt.(max.(lambda2_1, 0)))
    lambda1_1_smooth = smooth(lambda1_1)
    lambda2_1_smooth = smooth(lambda2_1)

    # Iterate until the relative change in f and g is small
    iterations = 0
    converged = false
    nNegative = 0
    while !converged && (iterations < functionIterateMethod.functionIterMax)
        iterations = iterations + 1

        fPrevious = fNew
        gPrevious = gNew

        # Estimate gradients (finite differences)
        fOld = smooth(fPrevious)
        gOld = smooth(gPrevious)
        fGrad = smooth(fdiffNU(xEvalPoints,fOld))
        gGrad = smooth(fdiffNU(xEvalPoints,gOld))

        # Find updated values
        fNew, gSquared = fixed_point_iterate(
            lambda1_1_smooth,
            lambda2_1_smooth,
            fOld,
            gOld,
            fGrad,
            gGrad,
            theta
        )
        nNegative += count(<(0), gSquared)
        gNew = sqrt.(max.(gSquared, 0))

        converged = relative_change(fNew, fPrevious) < tolerance &&
            relative_change(gNew, gPrevious) < tolerance
    end
    return fNew, gNew, iterations, converged, nNegative
end

# Change between iterates, relative to the size of the previous iterate
relative_change(new, previous) = sqrt(sum(abs2, new .- previous)/sum(abs2, previous))

# Finite difference iteration, giving updated f and g^2
function fixed_point_iterate(lambda1_1,lambda2_1,f,g,fGrad,gGrad,theta)
    fNew = lambda1_1 .- 0.5*g.*gGrad .- 
                0.5*theta*(fGrad.*g.*gGrad .- f.*gGrad.^2)
    gSquaredNew = lambda2_1 .- theta*(fGrad.*g.^2 .- f.*g.*gGrad)
    return fNew, gSquaredNew
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