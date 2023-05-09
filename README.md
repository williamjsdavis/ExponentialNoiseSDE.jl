# ExponentialNoiseSDE
Regression of stochastic processes with exponentially correleted noise. `ExponentialNoiseSDE` is a `julia` package designed to estimate a stochastic differential equation. The regression is flexable enough to account for white or exponentially correlated noise. It employs a block bootstrap algorithm to estimate uncertainties. 

# Notes

This is a port of the `MATLAB` package [`secn-PR`](https://github.com/williamjsdavis/secn-PR). See that repo for details of the theory and implementation.

An example of the use of this package is shown in the directory `example/`.

# Banchmarking

Individual steps in the method can be individually benchmarked. The result below shows that the main computational cost comes from calculation of the conditional moments.

<img src="/example/benchmarking_time.png" height="400"/>

# TODOs

- [ ] Add observations from multiple data sources
- [ ] Implement arbitrary `timeShiftSamplePoints` vector
- [ ] Implement arbitrary spatial evaluation points
- [x] Implement option of forcing a particular correlation time
- [ ] Add compatibility for combining observations of different time-steps
- [x] Make a plotting function for the distribution of correlation times
- [ ] Add more kernel functions to the library

# Changelog

- Version v0.1.0 - WIP

# References

- Lehle, B., & Peinke, J. (2018). Analyzing a stochastic process driven by Ornstein-Uhlenbeck noise. Physical Review E, 97(1), 012113.
- Lamouroux, D., & Lehnertz, K. (2009). Kernel-based regression of drift and diffusion coefficients of stochastic processes. Physics Letters A, 373(39), 3507-3512.
- Kunsch, H. R. (1989). The jackknife and the bootstrap for general stationary observations. The annals of Statistics, 1217-1241.