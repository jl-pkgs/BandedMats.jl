using Test
using LinearAlgebra

# speye(n) = SparseArrays.sparse(I, n, n)
ddmat(x::AbstractVector, d::Integer=2) = diff(diagm(x), d)

function WHIT(y::AbstractVector, w::AbstractVector; kw...)
  WHIT(y, w, 1:length(y); kw...)
end

function WHIT(y::AbstractVector, w::AbstractVector, x::AbstractVector;
  λ=2.0, p=2)
  D = ddmat(x, p)

  W = diagm(w)
  A = W + λ * D' * D

  L = cholesky(A).L # Matrix
  # L = cholesky(A, perm=1:n).L # sparse
  z = L' \ (L \ (w .* y))
  z
end


@testset "whit_band" begin
  n = 100
  y = rand(n)
  w = rand(n)
  x = rand(n)
  λ, p = 2.0, 3

  @time z1 = whit_band(y, w, x; λ=2.0, p=3)
  @time z2 = WHIT(y, w, x; λ=2.0, p=3)
  @test maximum(abs.(z1 - z2)) <= 1e-10
end
