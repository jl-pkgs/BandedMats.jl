@testset "LU_solve" begin
  # n, m = 8, 6 # 非齐次矩阵, LU分解
  n, m = 6, 6
  p, q = 2, 1
  A = rand(n, m)
  b = rand(n, 1)
  # force_band!(A, p, q)
  L, U = LU_full(A)

  r = A \ b
  r2 = solve_U(U, solve_L(L, b))
  @test r2 ≈ r

  r1 = U \ (L \ b)
  @test r1 ≈ r
end

@testset "LU_band_solve" begin
  n = 6
  p, q = 2, 1
  A = rand(n, n)
  b = rand(n, 1)
  force_band!(A, p, q)
  r = A \ b

  B = BandedMat(A, p, q; type="kong", zipped=false)
  BL, BU = LU_band(B)

  r2 = solve_U(BU, solve_L(BL, b))
  @test r2 ≈ r
end
