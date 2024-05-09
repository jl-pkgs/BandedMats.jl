using Test
using Symbolics
import Symbolics: scalarize, variables

Mat(name, n) = variables(name, 1:n, 1:n)
Vec(name, n) = variables(name, 1:n)

@testset "BandMatrix" begin
  p, q = 2, 1
  A = Mat(:a, 5)
  b = BandMatrix(A, p, q)
  bd1 = BandedMatrix(b; type="lapack")
  bd2 = BandedMatrix(b; type="kong")
  print(b)

  r1 = BandMatrix(bd1).data - A
  r2 = BandMatrix(bd2).data - A

  @test maximum(abs.(r1)) <= 1e-10
  @test maximum(abs.(r2)) <= 1e-10
end
