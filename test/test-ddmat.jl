@testset "ddmat" begin
  x = 1:10
  @test Matrix(ddmat_band(x, 5)) ≈ ddmat_full(x, 5)
  @test Matrix(ddmat_band(x, 3)) ≈ ddmat_full(x, 3)
  @test Matrix(ddmat_band(x, 1)) ≈ ddmat_full(x, 1)

  x = rand(10)
  @test Matrix(ddmat_band(x, 5)) ≈ ddmat_full(x, 5)
  @test Matrix(ddmat_band(x, 3)) ≈ ddmat_full(x, 3)
  @test Matrix(ddmat_band(x, 1)) ≈ ddmat_full(x, 1)
end
