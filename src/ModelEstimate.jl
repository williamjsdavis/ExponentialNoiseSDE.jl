# Estimating SDE model (drift and noise functions + correlation time)

# MATLAB equivalent: FitOptionsClass.m
struct ModelEstimateSettings
    thetaConvergenceValue::Float64
    functionConvergenceValue::Float64
    fixThetaFlag::Bool
    fixThetaValue::Float64
    keepObservation::Bool
    displayOutput::Bool
    # Default settings
    function ModelEstimateSettings(;
        thetaConvergenceValue = 1E-2,
        functionConvergenceValue = 0.2,
        fixThetaFlag = false,
        fixThetaValue = 0.0,
        keepObservation = true,
        displayOutput = true
    ) begin
        return new(
            thetaConvergenceValue,
            functionConvergenceValue,
            fixThetaFlag,
            fixThetaValue,
            keepObservation,
            displayOutput
        )
    end
    end
end


