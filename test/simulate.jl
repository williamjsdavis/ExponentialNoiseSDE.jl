## Ground-truth simulator for tests
#
# Simulates the model the package estimates,
#     dX/dt = f(X) + g(X) η(t),
#     dη    = -η/θ dt + (1/θ) dW,
# where η is Ornstein-Uhlenbeck noise with correlation time θ, normalised so
# that <η(t)η(s)> = exp(-|t-s|/θ)/(2θ). As θ → 0, η → white noise and X
# approaches the Stratonovich SDE dX = f dt + g ∘ dW.
#
# η is advanced with its exact Gaussian transition. For θ > 0 the X equation is
# an ODE driven by a continuous η, so X is advanced pathwise with Heun's method.
# Internal steps are no longer than θ/25.
#
# Paths can diverge: with the example noise g(x) = 1 + x²/4 a large excursion
# grows faster than the drift pulls it back. A diverged path throws, so callers
# never estimate from non-finite data.

# True drift and noise functions of the package example
example_drift(x) = -2.5*x
example_noise(x) = x > 0 ? 1 + x^2/4 : 1.0

function simulate_process(rng, f, g, θ, dt, N; burnin=2000)
    substeps = max(1, ceil(Int, 25*dt/θ))
    h = dt/substeps
    ρ = exp(-h/θ)
    σ = sqrt((1 - ρ^2)/(2θ))
    x = 0.0
    η = randn(rng)/sqrt(2θ) # Stationary start
    X = Vector{Float64}(undef, N)
    for n in (1 - burnin):N
        for _ in 1:substeps
            ηNew = ρ*η + σ*randn(rng)
            k1 = f(x) + g(x)*η
            xPredict = x + h*k1
            k2 = f(xPredict) + g(xPredict)*ηNew
            x += 0.5*h*(k1 + k2)
            η = ηNew
        end
        isfinite(x) || throw(DomainError(n, "simulated path diverged; choose another seed"))
        if n >= 1
            X[n] = x
        end
    end
    return X
end

# White-noise limit (θ = 0) with linear drift and unit noise, dX = -aX dt + dW,
# using the exact Ornstein-Uhlenbeck transition
function simulate_white_noise_process(rng, a, dt, N)
    ρ = exp(-a*dt)
    σ = sqrt((1 - ρ^2)/(2a))
    x = randn(rng)/sqrt(2a) # Stationary start
    X = Vector{Float64}(undef, N)
    for n in 1:N
        x = ρ*x + σ*randn(rng)
        X[n] = x
    end
    return X
end
