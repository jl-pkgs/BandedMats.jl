@testset "inv_diag" begin
  n = 10
  A = rand(n, n)
  p, q = 2, 2
  force_band!(A, p, q)
  force_sym!(A)

  # L * U = A
  l, d = LDL_full(A)
  bl = BandedMat(l, p, 0; zipped=false)
  @test inv_diag(bl, d) ≈ diag(inv(A))

  U = diagm(d) * l'
  U2 = BandedMat(U, 0, q; zipped=false).data[:, 2:end]
  @test inv_diag(U2, d; m=q) ≈ diag(inv(A))

  B = BandedL(A, p; zipped=false)
  BL, d = LDL_band(B)
  @test inv_diag(BL, d) ≈ diag(inv(A))
end
