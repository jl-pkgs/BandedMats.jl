using Test
using LinearAlgebra
using SparseArrays
using BandedMats

function Base.diff(x::SparseMatrixCSC, d::Integer=1)
  D = x[2:end, :] .- x[1:end-1, :]
  d >= 2 ? diff(D, d - 1) : D
end

# 相比于原版进行了微调
function ddmat(x::AbstractVector, d::Integer=2)
  n = length(x)
  if d == 0
    return sparse(I, n, n)
  else
    dx = (x[d+1:n] - x[1:(n-d)]) #./ d
    V = spdiagm(1 ./ dx)
    return V * diff(ddmat(x, d - 1))
  end
end

function WHIT(y::AbstractVector, w::AbstractVector; kw...)
  WHIT(y, w, 1:length(y); kw...)
end

function WHIT(y::AbstractVector, w::AbstractVector, x::AbstractVector;
  λ=2.0, p=2)
  n = length(y)
  D = ddmat(x, p)

  W = spdiagm(w)
  A = W + λ * D' * D

  # L = cholesky(A).L # Matrix
  L = cholesky(A, perm=1:n).L # sparse
  z = L' \ (L \ (w .* y))
  z
end
