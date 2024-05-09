# LDL矩阵分解算法
# https://en.wikipedia.org/wiki/Cholesky_decomposition

# SymBandedMat
function LDL(A::AbstractMatrix{T}, tol::Real=1e-10) where {T}
  n = size(A, 1)
  # L = SymTridiagonal(zeros(n), zeros(n - 1))
  L = zeros(T, n, n)
  D = zeros(T, n)

  for i = 1:n
    D[i] = A[i, i]
    for j = 1:i-1
      D[i] -= L[i, j]^2 * D[j]
    end
    abs(D[i]) < tol && error("LDL: matrix is not positive definite")

    L[i, i] = 1
    for j = i+1:n
      L[j, i] = A[j, i]
      for k = 1:i-1
        L[j, i] -= L[j, k] * L[i, k] * D[k]
      end
      L[j, i] /= D[i]
    end
  end
  return L, D
end

# 强制修改为对称
function force_sym!(A::AbstractMatrix{T}) where {T}
  n, m = size(A)
  @assert n == m "check_sym: matrix is not square"
  for i = 1:n
    for j = 1:i-1
      A[i, j] = A[j, 1]
    end
  end
  A
end

export LDL, force_sym!
