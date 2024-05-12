# using Pkg
# Pkg.activate(".")
include("main_whit.jl")

@testset "whit" begin
  n = 100
  y = rand(n)
  w = rand(n)
  # w = fill(1.0, n)

  x = 1:n
  λ, p = 2.0, 3

  # @time z1 = whit_band(y, w, x; λ=2.0, p=3)
  @time z = WHIT(y, w, x; λ=2.0, p=3)
  @time z3, cve = whit3(y, w; lambda=2.0, include_cve=false)
  @test maximum(abs.(z - z3)) <= 1e-10

  @time z = WHIT(y, w, x; λ=2.0, p=2)
  @time z2, cve2 = whit2(y, w; lambda=2.0, include_cve=false)
  @test maximum(abs.(z - z2)) <= 1e-10
end
