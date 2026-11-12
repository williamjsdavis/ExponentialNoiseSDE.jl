# Kernel Based Regression of conditional moments

# Sums over pairs (X[t], X[t+τ]) from which the conditional moments and the
# autocorrelation increments are computed. They add up over any set of pairs,
# so they can be accumulated in parts.
struct MomentSums
    weight::Matrix{Float64} # Σ K, time shifts × evaluation points
    increment::Matrix{Float64} # Σ K ΔX
    squaredIncrement::Matrix{Float64} # Σ K ΔX²
    unweightedSquaredIncrement::Vector{Float64} # Σ ΔX², one per time shift
    pairCount::Vector{Int64} # Number of pairs, one per time shift
end

function MomentSums(momentSettings::ConditionalMomentSettings)
    nτ = length(momentSettings.timeShiftSamplePoints)
    nx = length(momentSettings.evalPoints)
    return MomentSums(zeros(nτ,nx), zeros(nτ,nx), zeros(nτ,nx), zeros(nτ), zeros(Int64,nτ))
end

# Add the pairs (X[t], X[t+τ]) that start at t in starts, for each time shift
# τ in the settings with t+τ within the series
function add_pairs!(sums::MomentSums, X, starts, momentSettings::ConditionalMomentSettings)
    # Passing the kernel as an argument compiles the loop for its concrete type
    return add_pairs!(sums, X, starts, momentSettings.timeShiftSamplePoints, momentSettings.evalPoints,
        momentSettings.bandwidth, momentSettings.kernel)
end

function add_pairs!(sums::MomentSums, X, starts, timeShiftSamplePoints, xEvalPoints, bandwidth, kernel::Kernel)
    n = length(X)
    hinv = inv(bandwidth)
    kernel_scaled(x) = hinv*kernel(hinv*x)

    # Only evaluation points within the kernel support get nonzero weights.
    # The margin keeps points whose weight rounds to nonzero at the edge.
    supportHalfWidth = support_radius(kernel)*bandwidth*(1 + 1E-8)

    for i_left in starts
        X_left = X[i_left]
        for (i_ind,i_tau) in enumerate(timeShiftSamplePoints)
            i_right = i_left + i_tau
            i_right <= n || break # Shifts are increasing
            ΔX = X[i_right] - X_left
            sums.unweightedSquaredIncrement[i_ind] += ΔX*ΔX
            sums.pairCount[i_ind] += 1
        end
        j_first = searchsortedfirst(xEvalPoints, X_left - supportHalfWidth)
        j_last = searchsortedlast(xEvalPoints, X_left + supportHalfWidth)
        for j_ind in j_first:j_last
            # The weight depends on the sample and evaluation point, not the shift
            K_weight = kernel_scaled(xEvalPoints[j_ind] - X_left)
            for (i_ind,i_tau) in enumerate(timeShiftSamplePoints)
                i_right = i_left + i_tau
                i_right <= n || break
                ΔX = X[i_right] - X_left
                sums.weight[i_ind,j_ind] += K_weight
                sums.increment[i_ind,j_ind] += K_weight * ΔX
                sums.squaredIncrement[i_ind,j_ind] += K_weight * ΔX*ΔX
            end
        end
    end
    return sums
end

# Conditional mean and variance of the increments (NaN where the kernel weights
# are all zero), and the autocorrelation increments A(τ) - A(0). For a
# stationary process A(τ) - A(0) = -E[(X(t+τ) - X(t))²]/2, so they are minus
# half the mean squared increments.
function moments(sums::MomentSums)
    N = sums.weight
    #NOTE: No Bessel correction
    moment2 = (sums.squaredIncrement .- (sums.increment .* sums.increment) ./ N) ./ N
    moment1 = sums.increment ./ N
    dA = -0.5 .* sums.unweightedSquaredIncrement ./ sums.pairCount
    return moment1, moment2, dA
end
