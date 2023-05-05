# Sample autocorrelation

#NOTE: Single data only
function autocorr_increment(X,tauMax)
    acf = myAutocorr(X,tauMax)
    nX = length(X)
    varX = var(X)

    # Tidy up array
    dA = (acf[2:end]-acf[1])*varX

    return dA
end
function myAutocorr(X,lags)
    nX = length(X)
    Xdemean = X .- mean(X)
    nPower = nextpow(2, nX) + 1
    nFFT = 2^nPower
    Xpadded = zeros(nFFT)
    Xpadded[1:nX] = Xdemean
    F = fft(Xdemean)
    F = F .* conj(F)
    acf = ifft(F)
    acf = acf[1:(lags+1)]
    acf = acf ./ acf[1]
    return acf
end