## Bootstrap statistics for uncertainties
# MATLAB equivalent: SPbootstrapClass.m

struct BootstrapStatistics
    correlationEstimate
    driftEstimate
    noiseEstimate
    modelError
    bootstrapSettings::BootstrapSettings
end


# Pack bootstrap results into a dict
function generate_statistics_dict(varMean,varMedian,varStd,varPrc)
    return Dict(
        :mean=>varMean,
        :median=>varMedian,
        :std=>varStd,
        :percentiles95=>varPrc,
    )
end

# Calculate sample statistics from bootstrap samples
function calculate_sample_statistics(bootstrapEstimates::Array{BootstrapModelEstimate},bootstrapSettings::BootstrapSettings)
    # Samples
    corrSamples = broadcast(x->x.correlationEstimate, bootstrapEstimates)
    errorSamples = broadcast(x->x.modelError.meanAbsoluteError, bootstrapEstimates)
    get_drift_samples = j->broadcast(x->x.driftEstimate[j], bootstrapEstimates)
    get_noise_samples = j->broadcast(x->x.noiseEstimate[j], bootstrapEstimates)

    # Mean
    corrMean = mean(corrSamples)
    errorMean = mean(errorSamples)
    driftMean = broadcast(i->mean(get_drift_samples(i)),eachindex(bootstrapEstimates))
    noiseMean = broadcast(i->mean(get_noise_samples(i)),eachindex(bootstrapEstimates))

    # Median
    corrMedian = median(corrSamples)
    errorMedian = median(errorSamples)
    driftMedian = broadcast(i->mean(get_drift_samples(i)),eachindex(bootstrapEstimates))
    noiseMedian = broadcast(i->mean(get_noise_samples(i)),eachindex(bootstrapEstimates))
    
    # Standard deviations
    corrStd = std(corrSamples)
    errorStd = std(errorSamples)
    driftStd = broadcast(i->std(get_drift_samples(i)),eachindex(bootstrapEstimates))
    noiseStd = broadcast(i->std(get_noise_samples(i)),eachindex(bootstrapEstimates))

    # Percentiles
    setPrc = [0.025,0.975]
    corrPrc = quantile(corrSamples,setPrc)
    errorPrc = quantile(errorSamples,setPrc)
    driftPrc = broadcast(i->quantile(get_drift_samples(i),setPrc),eachindex(bootstrapEstimates))
    noisePrc = broadcast(i->quantile(get_noise_samples(i),setPrc),eachindex(bootstrapEstimates))

    return BootstrapStatistics(
        generate_statistics_dict(corrMean,corrMedian,corrStd,corrPrc),
        generate_statistics_dict(driftMean,driftMedian,driftStd,driftPrc),
        generate_statistics_dict(noiseMean,noiseMedian,noiseStd,noisePrc),
        generate_statistics_dict(errorMean,errorMedian,errorStd,errorPrc),
        bootstrapSettings
    )
end