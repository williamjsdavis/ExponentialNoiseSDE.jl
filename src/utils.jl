# Utility and statistic functions

# Moving average, inspired my MATLAB
"""
M = movmean(X,n) returns an array of local k-point mean values, where each mean is calculated over a sliding window of length k across neighboring elements of A. 
When k is odd, the window is centered about the element in the current position. 
When k is even, the window is centered about the current and previous elements. 
The window size is automatically truncated at the endpoints when there are not enough elements to fill the window. 
When the window is truncated, the average is taken over only the elements that fill the window. M is the same size as A.
"""
function movmean(x,k)
    N = length(x)
    out = zeros(N)
    if iseven(k)
        movmean_even!(out,x,k,N)
    else
        movmean_odd!(out,x,k,N)
    end
    return out
end

# When k is odd, the window is centered about the element in the current position.
function movmean_odd!(out,x,k,N)
    j = (k-1) ÷ 2
    for i in 1:N
        imin = max(1,i-j)
        imax = min(N,i+j)
        out[i] = sum(view(x,imin:imax))/(imax-imin+1)
    end
end

# When k is even, the window is centered about the current and previous elements. 
function movmean_even!(out,x,k,N)
    j = k ÷ 2
    for i in 1:N
        imin = max(1,i-j)
        imax = min(N,i+j-1)
        out[i] = sum(view(x,imin:imax))/(imax-imin+1)
    end
end
