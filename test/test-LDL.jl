using LinearAlgebra
using BandedMats
using Test

@testset "LDL" begin
  n = 10
  A = rand(n, n)
  p, q = 2, 2
  force_band!(A, p, q)
  force_sym!(A)
  b = rand(n)

  ## LDL_full
  L, d = LDL_full(A)
  mat_equal(L * diagm(d) * L', A)

  ## LDL_band
  # B = BandedMat(A, p, q; zipped=false)
  A2 = BandedL(A, p; zipped=false)
  BL, d2 = LDL_band(A2)

  @test BL.data ≈ band_zip(L, p, 0)[:, 1:end-1]
  @test d ≈ d2

  ## LCL_solve
  @test LDL_solve(BL, d, b) ≈ A \ b
end
