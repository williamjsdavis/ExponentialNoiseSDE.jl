# Estimating SDE model (drift and noise functions + correlation time)

# MATLAB equivalent: FitOptionsClass.m
"""
    ModelEstimateSettings(; thetaMethod=EstimateTheta(), functionIterateMethod=SmoothedFiniteDiff())

How to estimate the model: the correlation time ([`EstimateTheta`](@ref) or
[`FixedTheta`](@ref)), and the drift and noise solver ([`SmoothedFiniteDiff`](@ref)
or [`FiniteDiff`](@ref)).
"""
Base.@kwdef struct ModelEstimateSettings
    thetaMethod::ThetaMethod = EstimateTheta()
    functionIterateMethod::FunctionIterateMethod = SmoothedFiniteDiff()
end

# MATLAB equivalent: SPmodelClass.m
"""
    ModelEstimate

An estimated model, from [`estimate_model`](@ref). Use [`correlation_time`](@ref),
[`drift`](@ref), [`noise`](@ref) and [`eval_points`](@ref) for the estimates.
`functionDiagnostics` records the drift and noise iteration (iterations,
convergence, negative squared noise), and `modelError` the mean absolute error
of the fitted moments.
"""
struct ModelEstimate
    correlationEstimate::Float64
    driftEstimate::Array{Float64,1}
    noiseEstimate::Array{Float64,1}
    driftInitial::Array{Float64,1}
    noiseInitial::Array{Float64,1}
    functionDiagnostics::FunctionIterationDiagnostics
    modelError::ModelError
    conditionalMoments::ConditionalMoments
    thetaProperties::ThetaProperties
    lambdaProperties::LambdaProperties
    fitSettings::ModelEstimateSettings
end

"""
    correlation_time(modelEstimate::ModelEstimate) -> Float64
    correlation_time(bootstrapStatistics) -> BootstrapSummary

The correlation time θ of the noise, in units of time; for bootstrap
statistics, a summary (`mean`, `median`, `std`, `lower95`, `upper95`).
"""
correlation_time(modelEstimate::ModelEstimate) = modelEstimate.correlationEstimate

"""
    drift(modelEstimate::ModelEstimate) -> Vector{Float64}
    drift(bootstrapStatistics) -> BootstrapSummary

The drift f at the evaluation points ([`eval_points`](@ref)), in units of x per
unit time; for bootstrap statistics, a summary at each point.
"""
drift(modelEstimate::ModelEstimate) = modelEstimate.driftEstimate

"""
    noise(modelEstimate::ModelEstimate) -> Vector{Float64}
    noise(bootstrapStatistics) -> BootstrapSummary

The noise amplitude g at the evaluation points, in units of x per square root
of unit time; for bootstrap statistics, a summary at each point.
"""
noise(modelEstimate::ModelEstimate) = modelEstimate.noiseEstimate

"""
    eval_points(modelEstimate::ModelEstimate) -> Vector{Float64}

The evaluation points x at which the drift and noise are estimated.
"""
eval_points(modelEstimate::ModelEstimate) = modelEstimate.conditionalMoments.xEvalPoints

# MATLAB equivalent: estimateSPmodel.m
"""
    estimate_model(conditionalMoments::ConditionalMoments, fitSettings::ModelEstimateSettings) -> ModelEstimate
    estimate_model(observation::Observation, momentSettings, fitSettings) -> ModelEstimate

Estimate the correlation time θ, drift f and noise g from conditional moments
(or directly from an observation, building the moments first).
"""
function estimate_model(conditionalMoments::ConditionalMoments, fitSettings::ModelEstimateSettings)
    # Estimate correlation time (theta)
    correlationEstimate, thetaProperties = estimate_theta(conditionalMoments, fitSettings.thetaMethod)

    # Estimate drift and noise properties
    lambdaProperties = lambda_search_linear(
        conditionalMoments.moment1,
        conditionalMoments.moment2,
        thetaProperties.rMatrix
    )
    
    # Estimate drift and noise functions
    lambda1_1 = lambdaProperties.lambda1Est[1,:] # Array of lambda^(1)_1
    lambda2_1 = lambdaProperties.lambda2Est[1,:] # Array of lambda^(2)_1
    driftEstimate, noiseEstimate, driftInitial, noiseInitial, functionDiagnostics = fg_solve_with_data(
        lambda1_1,
        lambda2_1,
        correlationEstimate,
        conditionalMoments.xEvalPoints,
        fitSettings.functionIterateMethod
    )

    # Mean fit error
    modelError = mean_fit_error(
        conditionalMoments,
        thetaProperties,
        lambdaProperties
    )

    return ModelEstimate(
        correlationEstimate,
        driftEstimate,
        noiseEstimate,
        driftInitial,
        noiseInitial,
        functionDiagnostics,
        modelError,
        conditionalMoments,
        thetaProperties,
        lambdaProperties,
        fitSettings
    )
end

# Other call args
function estimate_model(observation,momentSettings,fitSettings)
    # Calculate moments
    conditionalMoments = build_moments(observation, momentSettings)
    
    # Estimate model
    modelEstimate = estimate_model(conditionalMoments, fitSettings)
    return modelEstimate
end

Base.show(io::IO, modelEstimate::ModelEstimate) = print(io,
    "ModelEstimate(θ = $(round(correlation_time(modelEstimate), sigdigits=4)), ",
    "$(length(eval_points(modelEstimate))) evaluation points)")

function Base.show(io::IO, ::MIME"text/plain", modelEstimate::ModelEstimate)
    x = eval_points(modelEstimate)
    diagnostics = modelEstimate.functionDiagnostics
    function finiteRange(v)
        finiteValues = filter(isfinite, v)
        isempty(finiteValues) && return (NaN, NaN)
        return round.(extrema(finiteValues), sigdigits=4)
    end
    println(io, "ModelEstimate")
    println(io, "  correlation time θ: ", round(correlation_time(modelEstimate), sigdigits=6))
    println(io, "  evaluation points:  ", length(x), " in [", first(x), ", ", last(x), "]")
    println(io, "  drift range:        ", finiteRange(drift(modelEstimate)))
    println(io, "  noise range:        ", finiteRange(noise(modelEstimate)))
    println(io, "  iteration:          ", diagnostics.converged ? "converged" : "not converged",
        " in ", diagnostics.iterations, " iteration(s)")
    print(io, "  mean fit error:     ", round(modelEstimate.modelError.meanAbsoluteError, sigdigits=4))
end
