# Time-series observations

# MATLAB equivalent: ObservationClass.m
#NOTE: Only do single observations first
struct Observation
    X::Array{Float64,1}
    dt::Float64
    N::Int64
end
Observation(X,dt) = Observation(X,dt,length(X))

