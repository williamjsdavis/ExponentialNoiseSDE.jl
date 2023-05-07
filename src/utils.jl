# Utility and statistic functions

#moving_average(vs,n) = [sum(@view vs[i:(i+n-1)])/n for i in 1:(length(vs)-(n-1))]

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
        @show (imin,imax)
        out[i] = sum(view(x,imin:imax))/(imax-imin+1)
    end
end

## To add to tests
@testset "Statistics and utility functions" begin
    @test movmean([1,2,3,4,5],3) .== [1.5,2.0,3.0,4.0,4.5]
    @test movmean([1,2,3,4,5],2) .== [1.0,1.5,2.5,3.5,4.5]
    @test movmean([1,2,3,4,5],4) .== [1.5,2.0,2.5,3.5,4.0]
    @test movmean([1,2,3,4],3) .== [1.5,2.0,3.0,3.5]
    @test movmean([10,2,3,2,4,8,21],4) .== [6.0,5.0,4.25,2.75,4.25,8.75,11.0]
    @test movmean([9,2,4,3,2,7,21],5) .== [5.0,4.5,4.0,3.6,4.7,8.25,10.0]
end