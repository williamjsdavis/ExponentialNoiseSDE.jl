# Test data

## `exampleData.txt`

50 000 samples of a simulated process. The file is identical to
`example/exampleDataB.txt` in the MATLAB package
[`secn-PR`](https://github.com/williamjsdavis/secn-PR), whose
`example/importExampleData.m` gives the true model:

- time step `dt = 0.005`
- correlation time `θ = 0.01`
- drift `f(x) = -2.5x`
- noise `g(x) = 1 + x²/4` for `x > 0`, and `g(x) = 1` otherwise

The script that simulated it is not in either repository.
`test/simulate.jl` simulates the same model.

## `*MATLAB.txt`

Reference outputs computed with MATLAB `secn-PR` from `exampleData.txt`, with
15 lags, 20 evaluation points on `[-1, 1]` and bandwidth 0.1. They were made
with a version of `secn-PR` that differs from the public one: the public
version has no smoothed drift and noise solver. Both solvers ran one
fixed-point iteration (`count_max = 1` in `fgIter`), so the tests compare with
`functionIterMax = 1`. `driftEstimateMATLAB.txt` and `noiseEstimateMATLAB.txt`
are from the smoothed solver, and the `NoSmooth` files from the unsmoothed one.
The script that wrote these files is not in either repository.
