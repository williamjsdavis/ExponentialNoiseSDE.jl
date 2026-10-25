## Helpers for scientific tests

# Run the full pipeline with the settings used throughout the scientific tests
function fit_model(X, dt;
        lags = 1:15,
        nEvalPoints = 20,
        evalLims = (-1.0, 1.0),
        bandwidth = 0.1,
        method = SmoothedFiniteDiff())
    momentSettings = ConditionalMomentSettings(
        collect(lags), nEvalPoints, evalLims, Epanechnikov(), bandwidth)
    fitSettings = ModelEstimateSettings(functionIterateMethod = method)
    return estimate_model(Observation(X, dt), momentSettings, fitSettings)
end

# Same method with a fixed number of iterations, so results don't depend on
# the stopping threshold
pinned(::SmoothedFiniteDiff, n) = SmoothedFiniteDiff(functionConvergenceValue = 0.0, functionIterMax = n)
pinned(::FiniteDiff, n) = FiniteDiff(functionConvergenceValue = 0.0, functionIterMax = n)

# Estimates from data with time measured in units c times smaller (e.g. s → ms),
# mapped back to the original units: θ scales by c, f by 1/c and g by 1/sqrt(c)
function time_rescaled(X, dt, method; c = 1000.0)
    m = fit_model(X, c*dt; method)
    return correlation_time(m)/c, drift(m) .* c, noise(m) .* sqrt(c)
end

# Estimates from data transformed x → a*x + b, mapped back to the original
# scale: θ is unchanged, f and g scale by a
function space_rescaled(X, dt, method; a = 10.0, b = 0.5)
    m = fit_model(a .* X .+ b, dt; method, evalLims = (b - a, b + a), bandwidth = 0.1*a)
    return correlation_time(m), drift(m) ./ a, noise(m) ./ a
end

# Does model m agree with the (θ, f, g) triple from a transformed fit?
function matches(m, (θ, f, g); rtol = 1e-6)
    return isapprox(θ, correlation_time(m); rtol) &&
        isapprox(f, drift(m); rtol) &&
        isapprox(g, noise(m); rtol)
end
