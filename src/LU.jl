export LU_gauss, LU_full, LU_band, LU_band_symmetry
export LU_band_full

# 高斯消元法
LU_gauss(A) = LU_gauss!(deepcopy(A))
function LU_gauss!(A::AbstractMatrix{T}) where {T}
  _n, _m = size(A)
  m = min(_n, _m)
  n = _n
  L = zeros(T, n, m) # 第二项不确定

  for i = 1:m
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
  (; L, U=A[1:m, :])
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
  _n, _m = size(A)
  m = min(_n, _m)
  n = _n
  l = zeros(T, n, m)
  u = zeros(T, m, _m)

  for i = 1:m
    l[i, i] = 1
    for j = i:_m
      u[i, j] = A[i, j] - sum(l[i, 1:i-1] .* u[1:i-1, j])
    end

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
