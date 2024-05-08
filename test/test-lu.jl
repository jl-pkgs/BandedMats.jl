using LinearAlgebra
using BandMatrices
using Test


@testset "LU_gauss" begin
  n = 4
  A = rand(n, n)
  l, u = lu(A, NoPivot())

  L, U = LU_gauss(A)
  @test maximum(abs.(L - l)) <= 1e-10
  @test maximum(abs.(U - u)) <= 1e-10

  L, U = LU_full(A)
  @test maximum(abs.(L - l)) <= 1e-10
  @test maximum(abs.(U - u)) <= 1e-10
end
