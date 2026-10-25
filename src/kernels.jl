# Kernels
#TODO: Add kernel scale into struct

abstract type Kernel end

# Kernels are zero for |x| beyond their support radius (Inf if unbounded)
support_radius(::Kernel) = Inf

# Epanechnikov, scaled to unit variance (support [-sqrt(5), sqrt(5)])
struct Epanechnikov <: Kernel
end

support_radius(::Epanechnikov) = sqrt(5.0)

#NOTE: Maybe change this to support [-1,1]?
const epan_scale = 3*sqrt(5)/100

function (k::Epanechnikov)(x)
    x2 = x*x
    if x2 < 5.0
        return epan_scale*(5.0 - x2)
    else
        return 0.0
    end
end
