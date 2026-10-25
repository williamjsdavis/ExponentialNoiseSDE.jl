## Bootstrap methods for uncertainties
# MATLAB equivalent: SPbootstrapClass.m

# MATLAB equivalent: BootstrapOptionsClass.m
Base.@kwdef struct BootstrapSettings
    blockLength::Int64 = 500
    nSamples::Int64 = 20
    biasCorrection::Bool = true
    function BootstrapSettings(blockLength, nSamples, biasCorrection)
        blockLength >= 1 || throw(ArgumentError(
            "blockLength must be at least 1, got $blockLength"))
        # Standard deviations need at least 2 samples
        nSamples >= 2 || throw(ArgumentError(
            "nSamples must be at least 2, got $nSamples"))
        return new(blockLength, nSamples, biasCorrection)
    end
end

# Model estimates from each bootstrap resample, one entry (or column) per resample
struct BootstrapSamples
    correlationEstimate::Vector{Float64}
    driftEstimate::Matrix{Float64} # Evaluation points × resamples
    noiseEstimate::Matrix{Float64} # Evaluation points × resamples
    meanAbsoluteError::Vector{Float64}
end

# Exported function
function estimate_bootstrap_statistics(modelEstimate::ModelEstimate,bootstrapSettings::BootstrapSettings;
        rng::AbstractRNG = default_rng())
    # A block as long as the series only rotates it, which underestimates
    # the uncertainty
    observation = modelEstimate.conditionalMoments.observation
    bootstrapSettings.blockLength < length(observation) || throw(ArgumentError(
        "blockLength ($(bootstrapSettings.blockLength)) must be shorter than the series ($(length(observation)) samples)"))

    @debug "Starting correlated process regression and bootstrap iterations" bootstrapSettings
    samples = bootstrap_samples(
        observation,
        modelEstimate.conditionalMoments.momentSettings,
        modelEstimate.fitSettings,
        bootstrapSettings,
        rng
    )
    return summarize_bootstrap(modelEstimate, samples, bootstrapSettings)
end

# Re-estimate the model, with the same settings, on block bootstrap resamples
# of the observation. Warnings from the fits are reported once per kind.
function bootstrap_samples(observation, momentSettings, fitSettings, bootstrapSettings, rng)
    nSamples = bootstrapSettings.nSamples
    nEvalPoints = momentSettings.nEvalPoints
    correlationEstimate = Vector{Float64}(undef, nSamples)
    driftEstimate = Matrix{Float64}(undef, nEvalPoints, nSamples)
    noiseEstimate = Matrix{Float64}(undef, nEvalPoints, nSamples)
    meanAbsoluteError = Vector{Float64}(undef, nSamples)

    # May be able to parallelize
    warningCounter = WarningCounter(current_logger())
    with_logger(warningCounter) do
        for i in 1:nSamples
            @debug "Bootstrap sample $i of $nSamples"
            indices = block_bootstrap_indices(rng, length(observation), bootstrapSettings.blockLength)
            estimate = estimate_model(
                Observation(observation.X[indices], observation.dt),
                momentSettings,
                fitSettings
            )
            correlationEstimate[i] = estimate.correlationEstimate
            driftEstimate[:,i] = estimate.driftEstimate
            noiseEstimate[:,i] = estimate.noiseEstimate
            meanAbsoluteError[i] = estimate.modelError.meanAbsoluteError
        end
    end
    for (id, count) in warningCounter.counts
        @warn "$count of $nSamples bootstrap fits: $(warningCounter.messages[id])"
    end

    return BootstrapSamples(correlationEstimate, driftEstimate, noiseEstimate, meanAbsoluteError)
end

# Indices of a circular block bootstrap resample of 1:N: blocks of blockLength
# consecutive indices from random starts, wrapping from N back to 1, joined and
# cut to length N
#NOTE: The joins between blocks create increments that are not in the series
function block_bootstrap_indices(rng, N, blockLength)
    blockStarts = rand(rng, 1:N, ceil(Int, N/blockLength))
    indices = [mod(start + k - 1, N) + 1 for start in blockStarts for k in 0:blockLength-1]
    return indices[1:N]
end

# Logger that counts warnings by the statement that logged them, so warnings
# repeated across bootstrap fits can be reported once. Other messages are
# passed to the parent logger.
struct WarningCounter <: AbstractLogger
    parent::AbstractLogger
    counts::Dict{Any,Int}
    messages::Dict{Any,String}
end
WarningCounter(parent) = WarningCounter(parent, Dict{Any,Int}(), Dict{Any,String}())

is_warning(level) = Logging.Warn <= level < Logging.Error

Logging.min_enabled_level(counter::WarningCounter) =
    min(Logging.Warn, Logging.min_enabled_level(counter.parent))
Logging.shouldlog(counter::WarningCounter, level, args...) =
    is_warning(level) || Logging.shouldlog(counter.parent, level, args...)
Logging.catch_exceptions(counter::WarningCounter) = Logging.catch_exceptions(counter.parent)

function Logging.handle_message(counter::WarningCounter, level, message, _module, group, id, file, line; kwargs...)
    if is_warning(level)
        counter.counts[id] = get(counter.counts, id, 0) + 1
        get!(counter.messages, id, string(message))
    else
        Logging.handle_message(counter.parent, level, message, _module, group, id, file, line; kwargs...)
    end
    return nothing
end
