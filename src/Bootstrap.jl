## Bootstrap methods for uncertainties
# MATLAB equivalent: SPbootstrapClass.m

# MATLAB equivalent: BootstrapOptionsClass.m
struct BootstrapSettings
    blockLength::Int64
    nSamples::Int64
    biasCorrection::Bool
    displayOutputFlag::Bool
end

# A single bootstrap model estimate
struct BootstrapModelEstimate
    correlationEstimate::Float64
    driftEstimate::Array{Float64,1}
    noiseEstimate::Array{Float64,1}
    modelError
end

# Exported function
function estimate_bootstrap_statistics(modelEstimate::ModelEstimate,bootstrapSettings::BootstrapSettings)
    # Toggle printing
    verbose(
        bootstrapSettings.displayOutputFlag, 
        "Starting correlated process regression and bootstrap iterations"
    )

    # Vary observation, keep the same settings
    estimate_model_set(obs) = 
        estimate_model(
            obs, 
            modelEstimate.conditionalMoments.momentSettings, 
            modelEstimate.fitSettings
        )
    
    originalObs = modelEstimate.conditionalMoments.obervation
    bootsrapEstimates = block_bootstrap(originalObs, estimate_model_set, bootstrapSettings)

    bootstrapStatistics = calculate_sample_statistics(modelEstimate,bootsrapEstimates,bootstrapSettings)

    return bootstrapStatistics
end

# Initialize results array and create bootstrap sampler
function block_bootstrap(obs,estimate_model_set,bootstrapSettings)
    bootstrapSampler = make_bootstrap_index_sampler(obs.N,bootstrapSettings.blockLength)
    bootsrapEstimates = Array{BootstrapModelEstimate}(undef, bootstrapSettings.nSamples)

    # May be able to parallelize
    dispatch_bootstrap_calculations!(
        bootsrapEstimates,
        obs,
        estimate_model_set,
        bootstrapSampler,
        bootstrapSettings.displayOutputFlag
    )
    return bootsrapEstimates
end

# Calculate all samples and estimates
function dispatch_bootstrap_calculations!(bootsrapEstimates,obs,estimate_model_set,bootstrapSampler,displayOutputFlag)
    for i in eachindex(bootsrapEstimates)
        verbose(
            displayOutputFlag, 
            "Sample: $i"
        )
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
function make_bootstrap_index_sampler(N,blockLength)
    rdraw2() = rand(1:N,ceil(Int, N/blockLength))
    drawin(i) = i:(i+blockLength-1)
    moddraw(i) = mod(i-1,N) + 1
    drawN(i) = moddraw.(i[1:N])
    allin(i) = reduce(vcat,broadcast(collect∘drawin,i))
    bsdrawi() = drawN(allin(rdraw2()))
    return bsdrawi
end



