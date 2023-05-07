## Bootstrap methods for uncertainties

#NOTE: Find out if this can be done in parallel
# MATLAB equivalent: BootstrapOptionsClass.m
struct BootstrapSettings
    blockLength::Int64
    nSamples::Int64
    displayOutputFlag::Bool
end

struct BootstrapModelEstimate
    correlationEstimate::Float64
    driftEstimate::Array{Float64,1}
    noiseEstimate::Array{Float64,1}
    modelError
end

# MATLAB equivalent: SPbootstrapClass.m
struct BootstrapEstimate
    distributions
    standardErrors
    percentiles95
    meanAbsoluteError
    modelEstimate::ModelEstimate
    bootstrapSettings::BootstrapSettings
end

# Estimate bootstrap uncertainties
function estimate_bootstrap_uncertainties(modelEstimate::ModelEstimate,bootstrapSettings::BootstrapSettings)
    # Code here
    verbose(message) = verbose(bootstrapSettings.displayOutputFlag,message)

    verbose("Starting correlated process regression and bootstrap iterations")
    #TODO: Add other messages

    # Vary observation, keep the same settings
    estimate_model_set(obs) = estimate_model(obs, momentSettings, fitSettings)
    

end

function block_bootstrap(obs,estimate_model_set,bootstrapSettings)
    bootstrapSampler = make_bootstrap_index_sampler(obs.N,bootstrapSettings.blockLength)

    for i in 1:bootstrapSettings.nSamples
        println(i)
        bootstrapModelEstimate = make_bootstrap_model(obs,bootstrapSampler)
    end
end

# Terrible anon. function-based implementation...
function make_bootstrap_index_sampler(N,blockLength)
    rdraw2() = rand(1:N,ceil(Int, N/blockLength))
    drawin(i) = i:(i+blockLength-1)
    moddraw(i) = mod(i-1,N) + 1
    drawN(i) = moddraw.(i[1:N])
    allin(i) = reduce(vcat,broadcast(collect∘drawin,i))
    bsdrawi() = drawN(allin(rdraw2()))
    return bsdrawi
end

#bootstrapObservation = make_bootstrap_sample(obs::Observation)
function make_bootstrap_sample(obs,bootstrapSampler)
    # Code here
    X_boostrap = obs.X[bootstrap_index_sampler()]
    return Observation(
        X_boostrap,
        obs.dt
    )
end

function make_bootstrap_model(obs,bootstrapSampler)
    bootstrapObs = make_bootstrap_sample(obs,bootstrapSampler)
    bootstrapModelEstimate = estimate_model_set(bootstrapObs)
    return BootstrapModelEstimate(
        bootstrapModelEstimate.correlationEstimate,
        bootstrapModelEstimate.driftEstimate,
        bootstrapModelEstimate.noiseEstimate,
        bootstrapModelEstimate.modelError
    )
end

