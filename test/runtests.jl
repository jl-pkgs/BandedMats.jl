using LinearAlgebra
using BandedMats
using Test

include("test-band.jl")
include("test-lu.jl")
include("test-solve.jl")


# # 代数余子式
# @testset "complement" begin
#   A = rand(4, 4)
#   R = complement(A)
#   @test size(R) == size(A)
# end
