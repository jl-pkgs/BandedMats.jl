@testset "inv_diag" begin
  n = 10
  A = rand(n, n)
  p, q = 2, 2
  force_band!(A, p, q)
  force_sym!(A)

  # U = diagm(d) * L'
  # L * U = A
  l, d = LDL_full(A)
  bl = BandedMat(l, p, 0; zipped=false)
  @test inv_diag(bl, d) ≈ diag(inv(A))

  B = BandedMat(A, p, q; zipped=false)
  BL, d = LDL_band(B)
  @test inv_diag(BL, d) ≈ diag(inv(A))
end
