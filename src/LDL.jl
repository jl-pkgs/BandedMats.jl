# LDL矩阵分解算法
# https://en.wikipedia.org/wiki/Cholesky_decomposition
"""
    LDL_band(B::BandedL{T}; tol::Real=1e-10)

# Example
```julia
n = 10
A = rand(n, n)
p, q = 2, 2
force_band!(A, p, q)
force_sym!(A)

B = BandedL(A, p; zipped=false)
BL, d = LDL_band(B)
```
"""
function LDL_band(B::BandedL{T}) where {T}
  n = size(B.data, 1)
  p = B.p

  L = zeros(T, n, p)
  BL = BandedL(L, p; size=(n, n))
  d = zeros(T, n)

  LDL_band!(BL, d, B)
end

function LDL_band!(BL::BandedL{T}, d::AbstractVector{T}, B::BandedL{T}) where {T}
  p::Int = B.p
  n::Int = size(B.data, 1)
  A::Matrix{T} = B.data
  L::Matrix{T} = BL.data

  # [i, j] => [i, j-i+p+1]
  @fastmath @inbounds for i = 1:n
    d[i] = A[i, p+1]
    for j = max(i - p, 1):i-1
      # 1 <= j-i+p+1 <= p
      d[i] -= L[i, j-i+p+1] * L[i, j-i+p+1] * d[j]
    end
    # abs(d[i]) < tol && error("LDL: matrix is not positive definite")
    for j = i+1:min(i + p, n)
      L[j, i-j+p+1] = A[j, i-j+p+1]

      for k = max(i - p, j - p, 1):i-1
        L[j, i-j+p+1] -= L[j, k-j+p+1] * L[i, k-i+p+1] * d[k]
      end
      L[j, i-j+p+1] /= d[i]
    end
  end
  return BL, d
end


# SymBandedMat
function LDL_full(A::AbstractMatrix{T}, tol::Real=1e-10) where {T}
  n, m = size(A)
  @assert n == m "LDL: matrix is not square"

  # L = SymTridiagonal(zeros(n), zeros(n - 1))
  L = zeros(T, n, n)
  d = zeros(T, n)

  @inbounds for i = 1:n
    d[i] = A[i, i]
    for j = 1:i-1
      d[i] -= L[i, j]^2 * d[j]
    end
    abs(d[i]) < tol && error("LDL: matrix is not positive definite")

    L[i, i] = 1
    for j = i+1:n
      L[j, i] = A[j, i]
      for k = 1:i-1
        L[j, i] -= L[j, k] * L[i, k] * d[k]
      end
      L[j, i] /= d[i]
    end
  end
  return L, d
end

export LDL_full, LDL_band
