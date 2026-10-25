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

# A single bootstrap model estimate
struct BootstrapModelEstimate
    correlationEstimate::Float64
    driftEstimate::Array{Float64,1}
    noiseEstimate::Array{Float64,1}
    modelError
end

# Exported function
function estimate_bootstrap_statistics(modelEstimate::ModelEstimate,bootstrapSettings::BootstrapSettings;
        rng::AbstractRNG = default_rng())
    # A block as long as the series only rotates it, which underestimates
    # the uncertainty
    originalObs = modelEstimate.conditionalMoments.obervation
    bootstrapSettings.blockLength < originalObs.N || throw(ArgumentError(
        "blockLength ($(bootstrapSettings.blockLength)) must be shorter than the series ($(originalObs.N) samples)"))

    @debug "Starting correlated process regression and bootstrap iterations" bootstrapSettings

    # Vary observation, keep the same settings
    estimate_model_set(obs) = 
        estimate_model(
            obs, 
            modelEstimate.conditionalMoments.momentSettings, 
            modelEstimate.fitSettings
        )

    bootsrapEstimates = block_bootstrap(originalObs, estimate_model_set, bootstrapSettings, rng)

    bootstrapStatistics = calculate_sample_statistics(modelEstimate,bootsrapEstimates,bootstrapSettings)

    return bootstrapStatistics
end

# Initialize results array and create bootstrap sampler
function block_bootstrap(obs,estimate_model_set,bootstrapSettings,rng)
    bootstrapSampler = make_bootstrap_index_sampler(obs.N,bootstrapSettings.blockLength,rng)
    bootsrapEstimates = Array{BootstrapModelEstimate}(undef, bootstrapSettings.nSamples)

    # May be able to parallelize
    dispatch_bootstrap_calculations!(
        bootsrapEstimates,
        obs,
        estimate_model_set,
        bootstrapSampler
    )
    return bootsrapEstimates
end

# Calculate all samples and estimates
function dispatch_bootstrap_calculations!(bootsrapEstimates,obs,estimate_model_set,bootstrapSampler)
    for i in eachindex(bootsrapEstimates)
        @debug "Bootstrap sample $i of $(length(bootsrapEstimates))"
        bootsrapEstimates[i] = make_bootstrap_estimate(obs,estimate_model_set,bootstrapSampler)
    end
end

# Create a single bootstrap sample and model estimate
function make_bootstrap_estimate(obs,estimate_model_set,bootstrapSampler)
    bootstrapObs = make_bootstrap_sample(obs,bootstrapSampler)
    bootstrapModelEstimate = estimate_model_set(bootstrapObs)
    return BootstrapModelEstimate(
        bootstrapModelEstimate.correlationEstimate,
        bootstrapModelEstimate.driftEstimate,
        bootstrapModelEstimate.noiseEstimate,
        bootstrapModelEstimate.modelError
    )
end

# Create a bootstrap observation
function make_bootstrap_sample(obs,bootstrapSampler)
    # Code here
    X_boostrap = obs.X[bootstrapSampler()]
    return Observation(
        X_boostrap,
        obs.dt
    )
end

# Create a bootstrap sampler
#NOTE: Terrible anon. function-based implementation...
function make_bootstrap_index_sampler(N,blockLength,rng)
    rdraw2() = rand(rng,1:N,ceil(Int, N/blockLength))
    drawin(i) = i:(i+blockLength-1)
    moddraw(i) = mod(i-1,N) + 1
    drawN(i) = moddraw.(i[1:N])
    allin(i) = reduce(vcat,broadcast(collect∘drawin,i))
    bsdrawi() = drawN(allin(rdraw2()))
    return bsdrawi
end



