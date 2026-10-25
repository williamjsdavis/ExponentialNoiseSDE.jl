# Changelog

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
