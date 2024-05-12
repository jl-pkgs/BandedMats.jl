using Test
using LinearAlgebra
using SparseArrays
using BandedMats

# function Base.diff(x::AbstractMatrix, d::Integer=1)
#   D = x[2:end, :] .- x[1:end-1, :]
#   d >= 2 ? diff(D, d - 1) : D
# end

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
    dx = (x[d+1:n] - x[1:(n-d)]) ./ d
    V = spdiagm(1 ./ dx)
    return V * diff(ddmat(x, d - 1))
  end
end

function WHIT(y::AbstractVector, w::AbstractVector; kw...)
  WHIT(y, w, 1:length(y); kw...)
end

function WHIT(y::AbstractVector, w::AbstractVector, x::AbstractVector;
  λ=2.0, p=2, include_cve=false)

  n = length(y)
  D = ddmat(x, p)

  W = spdiagm(w)
  A = W + λ * D' * D
  # A2 = Matrix(sparse(A))
  # L = cholesky(A2).L # Matrix
  L = cholesky(A, perm=1:n).L # sparse
  z = L' \ (L \ (w .* y))

  ## also include cve
  cve = -999.0
  if include_cve
    inv_L = Matrix(sparse(L))^-1 # 这一步可能消耗了较多的时间
    H = inv_L' * inv_L * W # 借力cholesky

    h = diag(H)
    r = @. (y - z) / (1 - h)
    cve = sqrt(sum(r .* r .* w) / sum(w))
  end
  z, cve
end
