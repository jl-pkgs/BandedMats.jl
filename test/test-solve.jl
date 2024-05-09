@testset "LU_solve" begin
  n = 6
  p, q = 2, 1
  A = rand(n, n)
  b = rand(n, 1)
  check_band!(A, p, q)

  r = A \ b
  L, U = LU_full(A)
  r1 = U \ (L \ b)
  r1 ≈ r

  r2 = solve_U(U, solve_L(L, b))
  r2 ≈ r
end

@testset "LU_band_solve" begin
  n = 6
  p, q = 2, 1
  A = rand(n, n)
  b = rand(n, 1)
  check_band!(A, p, q)
  r = A \ b

  L, U = LU_band(A; p, q)
  BL = BandedMatrix(L, p, 0)
  BU = BandedMatrix(U, 0, q)

  r2 = solve_U(BU, solve_L(BL, b))
  @test r2 ≈ r
end
