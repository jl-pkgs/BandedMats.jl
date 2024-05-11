using LinearAlgebra
using BandedMats
using Test

mat_equal(x, y) = @test maximum(abs.(x - y)) <= 1e-10

include("test-band.jl")
include("test-inv_diag.jl")
include("test-LDL.jl")
include("test-LU.jl")
include("test-solve.jl")

include("test-whit.jl")
# # 代数余子式
# @testset "complement" begin
#   A = rand(4, 4)
#   R = complement(A)
#   @test size(R) == size(A)
# end
