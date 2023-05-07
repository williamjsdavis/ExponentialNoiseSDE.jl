## Quantify errors of model fit

struct ModelError
    meanAbsoluteErrorMoment1
    meanAbsoluteErrorMoment2
    meanAbsoluteError
end

function mean_fit_error(
    conditionalMoments::ConditionalMoments,
    thetaProperties::ThetaProperties,
    lambdaProperties::LambdaProperties,
    displayOutputFlag)

    # Fit differences
    moment1Difference = conditionalMoments.moment1 - 
        thetaProperties.rMatrix * lambdaProperties.lambda1Est
    moment2Difference = conditionalMoments.moment2 - 
        thetaProperties.rMatrix * lambdaProperties.lambda2Est

    # Fit summed mean absolute errors
    meanAbsoluteErrorMoment1 = mean(abs.(moment1Difference))
    meanAbsoluteErrorMoment2 = mean(abs.(moment2Difference))
    meanAbsoluteError = 0.5*(meanAbsoluteErrorMoment1 + meanAbsoluteErrorMoment2)

    if displayOutputFlag
        #=
        println('Absolute fit error on moments')
        println(['Mean M^(1) fit error: ',
            sprintf('%0.4e',meanAbsoluteError.moment1)])
        println(['Mean M^(2) fit error: ',
            sprintf('%0.4e',meanAbsoluteError.moment2)])
        println(['Mean fit error: ',
            sprintf('%0.4e',meanAbsoluteError.bothMoments)])
        =#
    end
    
    return ModelError(
        meanAbsoluteErrorMoment1,
        meanAbsoluteErrorMoment2,
        meanAbsoluteError
    )
end
