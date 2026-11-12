# Changelog

## Unreleased

### Breaking changes

- `ConditionalMomentSettings` takes the evaluation points, `evalPoints`, instead of `nEvalPoints` and `evalLims`; any strictly increasing points can be used. For the previous evenly spaced grid, replace `nEvalPoints = 20, evalLims = (-1.0, 1.0)` with `evalPoints = range(-1.0, 1.0, length = 20)` (#43).
- `ConditionalMoments.observation` is now `observations`, a vector (#41).
- `BootstrapSettings.blockLength` counts windows of the series rather than samples (#38).

### Fixed

- The autocorrelation increments, to which θ is fitted, came from the sample autocorrelation, which divides the n − ν products at lag ν by n. This biased θ low for short series (−13% at 2 000 samples, −7% at 5 000). They are now minus half the mean squared increment, which is unbiased, and are computed in the same pass as the moments. θ on the example data changes by +0.2% (#36).
- The bootstrap joined blocks of the series end to end, so increments across the joins paired unrelated points. They acted like extra white noise: the resampled θ was 48% low at the default `blockLength = 500` (84% at 100), and the drift and noise summaries were biased too. Resamples are now made of blocks of the windows the estimate uses (Künsch, 1989), so every increment in a resample is an increment of the series. `blockLength` now counts windows, and must be less than the number of windows, the series length less the largest time shift (#38).

### Changed

- FFTW is no longer a dependency (#36).

### Performance

- Without the FFT, the moments of a series of 10⁷ samples take about half the time and allocate kilobytes instead of 2.8 GiB (#36).

### Added

- Uneven evaluation points, e.g. at quantiles of the data (#43).
- Several observations of the same process, with the same time step, can be combined: `build_moments` and `estimate_model` accept a vector of `Observation`s, pool the increments within each, and the bootstrap draws from the windows of all of them. Concatenating series instead biased θ low (−26% for 50 series of 1 000 samples) (#41).
- A benchmark script, `benchmark/scaling.jl`, that times the estimation against the series length and writes the figure in the README (#34).

## v0.2.0

Review-driven fixes and cleanup. Numbers in parentheses are GitHub issues.

### Breaking changes

- `ConditionalMomentSettings` takes a `Kernel` instance, e.g. `Epanechnikov()`, instead of a string; the string was ignored (#13). It also has a keyword constructor with `Epanechnikov()` as the default kernel (#23).
- `ModelEstimateSettings` has two settings: `thetaMethod` (`EstimateTheta()` or `FixedTheta(θ)`, replacing `fixThetaFlag`/`fixThetaValue`) and `functionIterateMethod`. `thetaConvergenceValue`, `keepObservation` and `displayOutputFlag` are removed; they had no effect (#23).
- `BootstrapSettings` no longer has `displayOutputFlag`; progress is logged at debug level (#23).
- `functionConvergenceValue` is a relative change between iterations (default `1E-4`), not an absolute sum of squares (default `0.2`), so default drift and noise estimates now come from a converged iteration (#17).
- Bootstrap statistics are typed summaries with fields `mean`, `median`, `std`, `lower95` and `upper95` instead of `Dict`s with `:percentiles95`, and include the resampled estimates in `samples` (#25).
- `ConditionalMoments.obervation` is renamed `observation`; `Observation` no longer has the field `N` (use `length(observation)`) (#27).
- Invalid inputs throw an `ArgumentError` (#11, #21).

### Fixed

- The autocorrelation was circular rather than linear (unpadded FFT), unlike MATLAB; results now match the MATLAB references to about 1e-6 or better (#5).
- The second correlation time search used `ceil(sqrt(θ)/dt)` lags, which depends on the units of time, crashed when the first search gave θ = 0, and gave unidentifiable fits for 3 or fewer lags. θ is now fitted once over all lags (#9).
- Time shift lags other than `1:ν_max` gave wrong drift and noise or errors; any increasing lags now work, and `estimate_model` no longer reads the raw data (#11).
- The exported `Epanechnikov` was undefined, and a broken stub method is removed (#13).
- Evaluation points without data no longer spread NaN to their neighbours, the fit error or the bootstrap (#15).
- The smoothed drift and noise iteration could not reach small tolerances, and convergence depended on the units; not converging is now reported, and negative squared noise is set to 0 with a warning instead of reflected with `abs` (#17).
- The bootstrap median of the drift and noise was their mean (#19).
- The MATLAB comparison of the unsmoothed solver compared different iteration counts (#7).

### Added

- Accessors `correlation_time`, `drift`, `noise` and `eval_points`, and compact printing of results (#27).
- `rng` keyword for `estimate_bootstrap_statistics` (#1).
- Drift and noise iteration diagnostics in `ModelEstimate.functionDiagnostics` (#17).
- Bootstrap warnings are reported once per kind, with a count (#25).
- Tests against data simulated from the model, with calibrated tolerances, including invariance under changes of units (#3); CI and `[compat]` bounds (#1).
- Docstrings, a documentation site, a runnable example and this changelog (#31).

### Performance

- Conditional moments are 6–8 times faster, with identical results (#29).

## v0.1.0

Initial port of the MATLAB package `secn-PR`.
