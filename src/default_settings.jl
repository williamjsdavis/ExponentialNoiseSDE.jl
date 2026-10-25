## Default estimation settings

# Model estimation
function ModelEstimateSettings(;
        thetaConvergenceValue = 1E-2,
        functionIterateMethod = SmoothedFiniteDiff(),
        fixThetaFlag = false,
        fixThetaValue = 0.0,
        keepObservation = true,
        displayOutputFlag = false
    ) begin
        return ModelEstimateSettings(
            thetaConvergenceValue,
            functionIterateMethod,
            fixThetaFlag,
            fixThetaValue,
            keepObservation,
            displayOutputFlag
        )
    end
end

# Finite difference function solution (convergence value is a relative change)
function FiniteDiff(;
        functionConvergenceValue = 1E-4,
        functionIterMax = 20
    ) begin FiniteDiff(
        functionConvergenceValue,
        functionIterMax
    )
    end
end

function SmoothedFiniteDiff(;
        nPointSmooth = 3,
        functionConvergenceValue = 1E-4,
        functionIterMax = 20
    ) begin SmoothedFiniteDiff(
        nPointSmooth,
        functionConvergenceValue,
        functionIterMax
    )
    end
end

# Bootstrapping
function BootstrapSettings(;
        blockLength = 500,
        nSamples = 20,
        biasCorrection = true,
        displayOutputFlag = true
    ) begin BootstrapSettings(
        blockLength,
        nSamples,
        biasCorrection,
        displayOutputFlag
    )
    end
end
