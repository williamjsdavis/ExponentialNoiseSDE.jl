# Sample autocorrelation

#NOTE: Single data only
function autocorr_increment(X,tauMax)
    acf = myAutocorr(X,tauMax)
    nX = length(X)
    varX = var(X)

    ####
    [acf,dA] = deal(cell(1,ndata)); % Preallocate
    [nX,varX] = deal(zeros(1,ndata));
    
    acf = myAutocorr(X,tauMax)
    nX(nd) = numel(X); % Number of points
    varX(nd) = var(X); % Variance
    dA = (acf(2:end)-acf(1))*var(X); % Tidy up array
    
    %% Merge data
    dAmat = cell2mat(dA); % Make into matrix
    dA_full = (dAmat*nX')/sum(nX);
    return dA_full
end
function myAutocorr(X,lags)
    Xdemean = X - mean(X);
    nFFT = 2^(nextpow2(length(Xdemean))+1);
    F = fft(Xdemean,nFFT);
    F = F.*conj(F);
    acf = ifft(F);
    acf = acf(1:(lags+1));
    acf = real(acf);
    acf = acf./acf(1);
    return acf
end