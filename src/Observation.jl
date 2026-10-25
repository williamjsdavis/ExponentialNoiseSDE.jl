# Time-series observations

# MATLAB equivalent: ObservationClass.m
#NOTE: Only do single observations first
struct Observation
    X::Array{Float64,1}
    dt::Float64
    N::Int64
    function Observation(X, dt, N)
        N == length(X) || throw(ArgumentError(
            "N = $N does not match the number of samples, $(length(X))"))
        (dt > 0 && isfinite(dt)) || throw(ArgumentError(
            "time step dt must be positive and finite, got $dt"))
        all(isfinite, X) || throw(ArgumentError(
            "X contains $(count(!isfinite, X)) non-finite value(s)"))
        return new(X, dt, N)
    end
end
Observation(X,dt) = Observation(X,dt,length(X))

