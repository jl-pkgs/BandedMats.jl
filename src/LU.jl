export LU_gauss, LU_full, LU_band, LU_band_symmetry
export LU_band_full

# 高斯消元法
LU_gauss(A) = LU_gauss!(deepcopy(A))
function LU_gauss!(A::AbstractMatrix{T}) where {T}
  n = size(A, 1)
  L = zeros(T, size(A))

  @inbounds for i = 1:n-1
    r1 = A[i, :]
    L[i, i] = 1
    # U[i, :] = r1
    for j = i+1:n
      f = A[j, i] / A[i, i]
      L[j, i] = f
      A[j, :] .= A[j, :] .- (f * r1)
      # L[:, 1] = A₁[j, :]
    end
  end
  L[n, n] = 1
  (; L, U=A)
end

"""
LU_full(A::AbstractArray)

Doolittle's method，可避免修改矩阵A。

# 对称矩阵
```math
A = L U 
A = L D L'
U = D L', L = U' D^-1
```
"""
function LU_full(A::AbstractMatrix{T}) where {T}
  n = size(A, 1)
  u = zeros(T, n, n)
  l = zeros(T, n, n)
  # u = tri_upper(:u, n)
  # l = tri_lower(:l, n)

  @inbounds for i = 1:n
    l[i, i] = 1
    for j = i:n
      u[i, j] = A[i, j] - sum(l[i, 1:i-1] .* u[1:i-1, j])
    end
    # 一般矩阵
    for i2 = i+1:n
      l[i2, i] = (A[i2, i] - sum(l[i2, 1:i-1] .* u[1:i-1, i])) / u[i, i]
    end
  end
  l, u
end

# ## 对称矩阵的福利: L = U' D^-1
# for i2 = i+1:min(i + p, i + q, n)
#   l[i2, i] = u[i, i2] / u[i, i]
# end
