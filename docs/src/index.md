# ExponentialNoiseSDE.jl

Estimate a stochastic differential equation driven by exponentially correlated
noise from a sampled time series, with bootstrap uncertainties. This is a Julia
port of the MATLAB package [`secn-PR`](https://github.com/williamjsdavis/secn-PR).

## The model

The observed process ``X(t)`` follows

```math
\frac{dX}{dt} = f(X) + g(X)\,\eta(t), \qquad d\eta = -\frac{\eta}{\theta}\,dt + \frac{1}{\theta}\,dW,
```

with drift ``f``, noise amplitude ``g``, and Ornstein-Uhlenbeck noise ``\eta`` of
correlation time ``\theta``. The noise is normalised so that
``\langle \eta(t)\eta(s) \rangle = e^{-|t-s|/\theta}/(2\theta)``. As
``\theta \to 0`` it becomes white noise, and the model becomes the Stratonovich
SDE ``dX = f\,dt + g \circ dW``.

Units: ``\theta`` is in units of time (those of the time step `dt`), ``f`` in
units of ``x`` per unit time, and ``g`` in units of ``x`` per square root of unit
time. Estimates transform accordingly when the units of time or ``x`` change.

## The method

Following Lehle & Peinke (2018):

1. [`build_moments`](@ref) estimates the conditional mean and variance of the
   increments ``\Delta X = X(t+\tau) - X(t)`` given ``X(t) = x``, by kernel
   regression at a grid of evaluation points ``x`` and time shifts ``\tau``. It
   also computes the autocorrelation increments of ``X`` at those shifts.
2. [`estimate_model`](@ref) fits ``\theta`` to the autocorrelation increments,
   expands the moments in basis functions ``r_k(\tau, \theta)`` by least
   squares, and solves for ``f`` and ``g`` by fixed-point iteration.
3. [`estimate_bootstrap_statistics`](@ref) repeats the estimate on circular
   block bootstrap resamples of the series.

## Quick start

The example data are 50 000 samples simulated with ``\mathrm{d}t = 0.005``,
``\theta = 0.01``, ``f(x) = -2.5x`` and ``g(x) = 1 + x^2/4`` for ``x > 0``
(``g = 1`` otherwise).

```@example quickstart
using ExponentialNoiseSDE, DelimitedFiles, Random

X = vec(readdlm(joinpath(pkgdir(ExponentialNoiseSDE), "test", "data", "exampleData.txt")))
observation = Observation(X, 0.005)

momentSettings = ConditionalMomentSettings(
    timeShiftSamplePoints = 1:15,
    nEvalPoints = 20,
    evalLims = (-1.0, 1.0),
    bandwidth = 0.1,
)
modelEstimate = estimate_model(observation, momentSettings, ModelEstimateSettings())
```

```@example quickstart
correlation_time(modelEstimate)
```

```@example quickstart
[eval_points(modelEstimate) drift(modelEstimate) noise(modelEstimate)]
```

Uncertainties from the bootstrap, with a fixed random number generator for
reproducibility:

```@example quickstart
bootstrapStatistics = estimate_bootstrap_statistics(
    modelEstimate, BootstrapSettings(nSamples = 50); rng = Xoshiro(1))
```

```@example quickstart
driftSummary = drift(bootstrapStatistics)
[eval_points(modelEstimate) driftSummary.lower95 driftSummary.upper95]
```

The resampled estimates are in `bootstrapStatistics.samples`.

## Options

- **Correlation time:** `ModelEstimateSettings(thetaMethod = FixedTheta(θ))`
  uses a known ``\theta`` instead of estimating it.
- **Drift and noise solver:** [`SmoothedFiniteDiff`](@ref) (the default) or
  [`FiniteDiff`](@ref), with a relative convergence tolerance. Iterations and
  convergence are recorded in `modelEstimate.functionDiagnostics`.
- **Progress messages:** the bootstrap logs progress at debug level; enable it
  with `ENV["JULIA_DEBUG"] = "ExponentialNoiseSDE"`.

## Known limitations

- **Bootstrap block joins.** Resamples join blocks of the series end to end,
  so increments that span a join are not increments of the original process.
  They act like extra white noise and bias the resampled estimates. The
  resampled ``\theta`` is lower than the estimate from the full series,
  increasingly so for short blocks (about −50% at `blockLength = 500` on the
  example data). With `biasCorrection = true` the 95% intervals are shifted to
  centre on the original estimate, but their width and the reported means and
  medians still include this effect.
- **Evaluation points without data.** Where no sample lies within the kernel
  support of an evaluation point, the moments, drift and noise are NaN (with a
  warning); the other points are unaffected.
- **Negative squared noise.** If the fixed-point iteration gives ``g^2 < 0`` at
  a point, it is set to 0 with a warning. This happens mainly for short or
  sparse series.
- **Kernel bias.** Kernel regression is biased by about ``h^2`` times the
  curvature of the moments and the density gradient, which does not shrink
  with more data. Choose the bandwidth ``h`` with this in mind.
- **Single series, constant time step.** Observations from several sources or
  with different time steps are not yet supported.

## References

- Lehle, B., & Peinke, J. (2018). Analyzing a stochastic process driven by
  Ornstein-Uhlenbeck noise. *Physical Review E*, 97(1), 012113.
- Lamouroux, D., & Lehnertz, K. (2009). Kernel-based regression of drift and
  diffusion coefficients of stochastic processes. *Physics Letters A*, 373(39),
  3507–3512.
- Künsch, H. R. (1989). The jackknife and the bootstrap for general stationary
  observations. *The Annals of Statistics*, 1217–1241.
