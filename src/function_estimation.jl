# Estimation of functions f and g

# Function estimation methods
abstract type FunctionIterateMethod end

struct FiniteDiff <: FunctionIterateMethod
    functionConvergenceValue::Float64
    functionIterMax::Int64
end

struct SmoothedFiniteDiff 
    nPointSmooth::Int64
    functionConvergenceValue::Float64
    functionIterMax::Int64
end

# Finite difference solution for functions
function fg_solve(lambda1_1,lambda2_1,theta,xEvalPoints,functionIterateMethod::FiniteDiff)
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

function fg_solve(lambda1_1,lambda2_1,theta,xEvalPoints,functionIterateMethod::SmoothedFiniteDiff)
    # Starting values
    fInitial = lambda1_1
    gInitial = sqrt.(abs.(lambda2_1))
    fNew = movmean(fInitial,functionIterateMethod.nPointSmooth)
    gNew = movmean(gInitial,functionIterateMethod.nPointSmooth)
    lambda1_1_smooth = movmean(lambda1_1,functionIterateMethod.nPointSmooth)
    lambda2_1_smooth = movmean(lambda2_1,functionIterateMethod.nPointSmooth)

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
            lambda1_1_smooth,
            lambda2_1_smooth,
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

# Finite difference iteration
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