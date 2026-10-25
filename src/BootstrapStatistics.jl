## Bootstrap statistics for uncertainties
# MATLAB equivalent: SPbootstrapClass.m

# Summary of the bootstrap samples of an estimate. With bias correction, the
# 95% interval comes from the samples shifted so their mean is the original
# estimate; the mean, median and standard deviation are of the samples as drawn.
struct BootstrapSummary{T}
    mean::T
    median::T
    std::T
    lower95::T
    upper95::T
end

struct BootstrapStatistics
    correlationEstimate::BootstrapSummary{Float64}
    driftEstimate::BootstrapSummary{Vector{Float64}}
    noiseEstimate::BootstrapSummary{Vector{Float64}}
    modelError::BootstrapSummary{Float64}
    samples::BootstrapSamples
    bootstrapSettings::BootstrapSettings
end

# Summarize the samples of each estimate
function summarize_bootstrap(
        modelEstimate::ModelEstimate,
        samples::BootstrapSamples,
        bootstrapSettings::BootstrapSettings)
    biasCorrection = bootstrapSettings.biasCorrection
    return BootstrapStatistics(
        summarize(samples.correlationEstimate, modelEstimate.correlationEstimate, biasCorrection),
        summarize(samples.driftEstimate, modelEstimate.driftEstimate, biasCorrection),
        summarize(samples.noiseEstimate, modelEstimate.noiseEstimate, biasCorrection),
        summarize(samples.meanAbsoluteError, NaN, false),
        samples,
        bootstrapSettings
    )
end

# Summary of the samples of one value. Uses the finite samples, since
# resamples with no data at an evaluation point give NaN there.
function summarize(samples::AbstractVector, original::Real, biasCorrection::Bool)
    finiteSamples = filter(isfinite, samples)
    isempty(finiteSamples) && return BootstrapSummary(NaN, NaN, NaN, NaN, NaN)

    sampleMean = mean(finiteSamples)
    # Bias correction shifts the samples so their mean is the original estimate
    intervalSamples = biasCorrection ? finiteSamples .- (sampleMean - original) : finiteSamples
    lower95, upper95 = all(isfinite, intervalSamples) ?
        quantile(intervalSamples, [0.025, 0.975]) : (NaN, NaN)

    return BootstrapSummary(sampleMean, median(finiteSamples), std(finiteSamples), lower95, upper95)
end

# Summary at each evaluation point, from samples with one row per point
function summarize(samples::AbstractMatrix, original::AbstractVector, biasCorrection::Bool)
    pointSummaries = [summarize(samples[j,:], original[j], biasCorrection) for j in axes(samples, 1)]
    return BootstrapSummary(
        [s.mean for s in pointSummaries],
        [s.median for s in pointSummaries],
        [s.std for s in pointSummaries],
        [s.lower95 for s in pointSummaries],
        [s.upper95 for s in pointSummaries]
    )
end
