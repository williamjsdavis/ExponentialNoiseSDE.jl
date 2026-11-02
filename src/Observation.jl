# Time-series observations

# MATLAB equivalent: ObservationClass.m
#NOTE: Only do single observations first
"""
    Observation(X, dt)

A time series `X` sampled at a constant time step `dt`. `X` must be finite and
`dt` positive. `length(observation)` is the number of samples.
"""
struct Observation
    X::Array{Float64,1}
    dt::Float64
    function Observation(X, dt)
        (dt > 0 && isfinite(dt)) || throw(ArgumentError(
            "time step dt must be positive and finite, got $dt"))
        all(isfinite, X) || throw(ArgumentError(
            "X contains $(count(!isfinite, X)) non-finite value(s)"))
        return new(X, dt)
    end
end

# Number of samples
Base.length(observation::Observation) = length(observation.X)

Base.show(io::IO, observation::Observation) =
    print(io, "Observation($(length(observation)) samples, dt = $(observation.dt))")
