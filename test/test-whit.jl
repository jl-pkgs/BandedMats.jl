# using Pkg
# Pkg.activate(".")
include("main_whit.jl")

@testset "whit" begin
  n = 100
  y = rand(n)
  w = rand(n)
  # w = fill(1.0, n)

  x = 1.0:n
  λ, p = 2.0, 3

  # whit_band
  @time z, cve = WHIT(y, w, x; λ=2.0, p=3, include_cve=true)
  @time z_band = whit_band(y, w, x; λ=2.0, p=3)
  @time z_band2 = whit_band(y, w; λ=2.0, p=3)
  @test z ≈ z_band
  @test z ≈ z_band2

  # whit3
  @time z3, cve3 = whit3(y, w; λ=2.0, include_cve=true)
  @test maximum(abs.(z - z3)) <= 1e-10
  @test cve3 ≈ cve

  # whit2
  @time z, cve = WHIT(y, w, x; λ=2.0, p=2)
  @time z2, cve2 = whit2(y, w; λ=2.0, include_cve=false)
  @test maximum(abs.(z - z2)) <= 1e-10
  @test cve ≈ cve2
end
