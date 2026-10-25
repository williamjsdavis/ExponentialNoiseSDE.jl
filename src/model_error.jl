## Quantify errors of model fit

struct ModelError
    meanAbsoluteErrorMoment1
    meanAbsoluteErrorMoment2
    meanAbsoluteError
end

function mean_fit_error(
    conditionalMoments::ConditionalMoments,
    thetaProperties::ThetaProperties,
    lambdaProperties::LambdaProperties)

    # Fit differences
    moment1Difference = conditionalMoments.moment1 - 
        thetaProperties.rMatrix * lambdaProperties.lambda1Est
    moment2Difference = conditionalMoments.moment2 - 
        thetaProperties.rMatrix * lambdaProperties.lambda2Est

    # Fit summed mean absolute errors, over the evaluation points with data
    # (moments are NaN where there are none)
    meanAbsoluteErrorMoment1 = mean(abs, filter(isfinite, moment1Difference))
    meanAbsoluteErrorMoment2 = mean(abs, filter(isfinite, moment2Difference))
    meanAbsoluteError = 0.5*(meanAbsoluteErrorMoment1 + meanAbsoluteErrorMoment2)

    return ModelError(
        meanAbsoluteErrorMoment1,
        meanAbsoluteErrorMoment2,
        meanAbsoluteError
    )
end
