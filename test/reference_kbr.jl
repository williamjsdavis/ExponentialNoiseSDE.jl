## Naive reference implementation of the conditional moments
# The package's implementation before it skipped zero kernel weights and
# computed each weight once per sample (kept to check it gives the same moments)
function reference_kbr_moments(X, timeShiftSamplePoints, xEvalPoints, bandwidth, kernel)
    n = length(X)
    nτ = length(timeShiftSamplePoints)
    nx = length(xEvalPoints)
    N = zeros(nτ,nx)
    M1 = zeros(nτ,nx)
    M2 = zeros(nτ,nx)

    hinv = inv(bandwidth)
    kernel_scaled(x) = hinv*kernel(hinv*x)

    for (i_left, X_left) in enumerate(X[1:end-1])
        for (i_ind,i_tau) in enumerate(timeShiftSamplePoints)
            i_right = i_left + i_tau
            if i_right <= n
                ΔX = X[i_right] - X_left
                for (j_ind,x_eval) in enumerate(xEvalPoints)
                    K_weight = kernel_scaled(x_eval - X_left)
                    N[i_ind,j_ind] += K_weight
                    M1[i_ind,j_ind] += K_weight * ΔX
                    M2[i_ind,j_ind] += K_weight * ΔX*ΔX
                end
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