using Test
# using Symbolics
# import Symbolics: scalarize, variables

# Mat(name, n) = variables(name, 1:n, 1:n)
# Vec(name, n) = variables(name, 1:n)

@testset "BandMat" begin
  p, q = 2, 1
  # A = Mat(:a, 5)
  A = rand(5, 4)
  b = BandMat(A, p, q)
  bd1 = BandedMat(b; type="lapack")
  bd2 = BandedMat(b; type="kong")
  print(b)

  @test band_zip(b; type="kong") == bd2.data

  r1 = BandMat(bd1).data - A
  r2 = BandMat(bd2).data - A

  @test maximum(abs.(r1)) <= 1e-10
  @test maximum(abs.(r2)) <= 1e-10
end

@testset "force" begin
  A = rand(5, 5)
  force_lower!(A) # 注意会修改A
  force_upper!(A) 
  @test A == diagm(diag(A))
end

@testset "transpose" begin
  A = rand(5, 5)
  p, q = 1, 2
  b = BandedMat(A, p, q; zipped=false)
  bt = BandedMat(A', q, p; zipped=false)
  mat_equal((b').data, bt.data)

  A = rand(5, 5)
  p, q = 2, 1
  b = BandedMat(A, p, q; zipped=false)
  bt = BandedMat(A', q, p; zipped=false)
  mat_equal((b').data, bt.data)
end

@testset "mult" begin
  function test_mult(x, y)
    x2 = BandedMat(x)
    y2 = BandedMat(y)
    @test x * y ≈ x.data * y.data
    @test x2 * y2 ≈ x.data * y.data
  end

  x = BandMat(rand(5, 4), 1, 2)
  y = BandMat(rand(4, 10), 2, 3)
  test_mult(x, y)

  x = BandMat(rand(5, 4), 1, 2)
  y = BandMat(rand(4, 5), 2, 3)
  test_mult(x, y)
end
