## Estimation of lambda and other properties

# Lambda properties
struct LambdaProperties
    lambda1Est # Array{Float64,2} ?
    lambda2Est # Array{Float64,2} ?
end

# Lambda estimation
function lambda_search_linear(moment1,moment2,rMatrix)
    lambda1 = rMatrix \ moment1
    lambda2 = rMatrix \ moment2
    return LambdaProperties(lambda1, lambda2)
end
