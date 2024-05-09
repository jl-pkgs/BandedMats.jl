using LinearAlgebra
using BandedMats
using Test

mat_equal(x, y) = @test maximum(abs.(x - y)) <= 1e-10

@testset "LDL" begin
  n = 10
  A = rand(10, 10)
  p, q = 2, 2
  force_band!(A, p, q)
  force_sym!(A)
  A

  L, d = LDL_full(A)
  mat_equal(L * diagm(d) * L', A)

  B = BandedMat(A, p, q; zipped=false)
  BL, d2 = LDL_band(B)

  @test BL.data ≈ band_zip(L, p, 0)[:, 1:end-1]
  @test d ≈ d2
end
