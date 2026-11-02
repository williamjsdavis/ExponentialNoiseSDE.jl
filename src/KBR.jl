# Kernel Based Regression of conditional moments

function calculate_xEvalPoints(evalLims::Tuple{Float64,Float64}, nEvalPoints::Int64)
    return range(evalLims[1], evalLims[2], length = nEvalPoints)
end

# Kernel Based Regression (kbr) for conditional moments
function kbr_moments(X, timeShiftSamplePoints, xEvalPoints, bandwidth, kernel::Kernel)
    n = length(X)
    nτ = length(timeShiftSamplePoints)
    nx = length(xEvalPoints)
    N = zeros(nτ,nx)
    M1 = zeros(nτ,nx)
    M2 = zeros(nτ,nx)

    hinv = inv(bandwidth)
    kernel_scaled(x) = hinv*kernel(hinv*x)

    # Only evaluation points within the kernel support get nonzero weights.
    # The margin keeps points whose weight rounds to nonzero at the edge.
    supportHalfWidth = support_radius(kernel)*bandwidth*(1 + 1E-8)

    for i_left in 1:n-1
        X_left = X[i_left]
        j_first = searchsortedfirst(xEvalPoints, X_left - supportHalfWidth)
        j_last = searchsortedlast(xEvalPoints, X_left + supportHalfWidth)
        for j_ind in j_first:j_last
            # The weight depends on the sample and evaluation point, not the shift
            K_weight = kernel_scaled(xEvalPoints[j_ind] - X_left)
            for (i_ind,i_tau) in enumerate(timeShiftSamplePoints)
                i_right = i_left + i_tau
                i_right <= n || break # Shifts are increasing
                ΔX = X[i_right] - X_left
                N[i_ind,j_ind] += K_weight
                M1[i_ind,j_ind] += K_weight * ΔX
                M2[i_ind,j_ind] += K_weight * ΔX*ΔX
            end
        end
    end

    #NOTE: No Bessel correction
    for i in eachindex(N)
        M2[i] = (M2[i] - (M1[i]*M1[i])/N[i]) / N[i]
        M1[i] = M1[i] / N[i]
    end

    return M1, M2
end