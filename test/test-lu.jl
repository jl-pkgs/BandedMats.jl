using LinearAlgebra
using BandMatrices
using Test

mat_equal(x, y) = maximum(abs.(x - y)) <= 1e-10

@testset "LU_gauss" begin
  n = 4
  A = rand(n, n)
  l, u = lu(A, NoPivot())

  L, U = LU_gauss(A)
  @test mat_equal(L, l)
  @test mat_equal(U, u)
  
  L, U = LU_full(A)
  @test mat_equal(L, l)
  @test mat_equal(U, u)
end

@testset "LU_band" begin
  n = 6
  p, q = 2, 1
  A = rand(n, n)
  check_band!(A, p, q)
  
  # 标准答案
  l, u = lu(A, NoPivot())
  u2 = band_zip(u, 0, q; type="kong")
  l2 = band_zip(l, p, 0; type="kong")[:, 1:end-1]

  # 我的版本
  L, U = LU_band(A; p, q)
  @test U ≈ u2
  @test L ≈ l2

  _l, _u = LU_band_full(A; p, q)
  @test _l ≈ l
  @test _u ≈ u
end
