using ExponentialNoiseSDE
using Test

@testset "First tests" begin
    @test 1==1
end

@testset "Observations" begin
    obs = Observation([0.0,0.0,0.0],0.1)
    @test obs.dt == 0.1
end