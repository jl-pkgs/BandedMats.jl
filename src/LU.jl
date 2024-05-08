using LinearAlgebra
# SparseArrays

# 高斯消元法
function LU_gauss(A)
  n = size(A, 1)
  T = typeof(A)
  L = T(diagm(ones(n)))

  for i = 1:n-1
    r1 = A[i, :]
    # U[i, :] = r1
    for j = i+1:n
      f = A[j, i] / A[i, i]
      L[j, i] = f
      A[j, :] .= A[j, :] .- (f * r1)
      # 为啥要引入U，这样已经求解完成了
      # L[:, 1] = A₁[j, :]
    end
  end
  (; L, U=A)
end

# Doolittle's method
# > 可避免修改矩阵A
function LU_full(a; symmetry=false)
  n = size(a, 1)
  u = tri_upper(:u, n)
  l = tri_lower(:l, n)

  for i = 1:n
    l[i, i] = 1
    for j = i:n
      u[i, j] = a[i, j] - sum(l[i, 1:i-1] .* u[1:i-1, j])
    end

    if symmetry
      ## 对称矩阵的福利: L = U' D^-1
      # A = L U 
      # A = L D L'
      # U = D L', L = U' D^-1
      for i2 = i+1:min(i + p, i + q, n)
        l[i2, i] = u[i, i2] / u[i, i]
      end
    else
      # 一般矩阵
      for i2 = i+1:n
        l[i2, i] = (a[i2, i] - sum(l[i2, 1:i-1] .* u[1:i-1, i])) / u[i, i]
      end
    end
  end
  l, u
end

# 带状矩阵的LU分解, Doolittle's method
# - [x] 有效的减少循环次数
# - [x] 压缩数据存储
function LU_band(a; p=2, q=1, symmetry=false)
  p != q && (symmetry = false)
  n = size(a, 1)
  # u = tri_upper(:u, n)
  # l = tri_lower(:l, n)
  l = variables(:l, 1:p, 1:n)   # [i+k, i] -> [k, i]  ;  [i, j] -> [i-j, j]
  u = variables(:l, 1:n, 1:q+1) # [i, i+k] -> [i, k+1];  [i, j] -> [i, j-i+1]

  fill!(u, 0)
  fill!(l, 0)
  for i = 1:n
    ## U矩阵同时压缩数据
    for j = i:min(i + q, n)
      u[i, j-i+1] = a[i, j]
      for k = max(i - p, j - q, 1):min(i - 1, j - 1)
        # 1 <= j - k <= q  ==> j - q <= k <= j - 1
        # 1 <= i - k <= p  ==> i - p <= k <= i - 1
        u[i, j-i+1] -= l[i-k, k] * u[k, j-k+1]
      end
    end

    if symmetry
      # 对称矩阵的福利: L = U' D^-1
      for i2 = i+1:min(i + p, i + q, n)
        # 0 <= i2 - i <= q ==> i <= i2 <= q+i
        l[i2-i, i] = u[i, i2-i+1] / u[i, 1] # u[i, i2]
      end

    else
      # 一般带状矩阵
      for i2 = i+1:min(i + p, n)
        l[i2-i, i] = a[i2, i]
        # 1 <= i - k <= q  ==> i - q <= k <= i - 1
        # 1 <= i2 - k <= p  ==> i2 - p <= k <= i2 - 1
        for k = max(i2 - p, i - q, 1):min(i - 1, i2 - 1)
          l[i2-i, i] -= l[i2-k, k] * u[k, i-k+1]
        end
        l[i2-i, i] /= u[i, 1]
      end

    end # endif symmetry
  end
  l, u
end


function LU_band_symmetry(a; p=2, q=1)
  n = size(a, 1)
  l = variables(:l, 1:p, 1:n)   # [i+k, i] -> [k, i]  ;  [i, j] -> [i-j, j]
  u = variables(:l, 1:n, 1:q+1) # [i, i+k] -> [i, k+1];  [i, j] -> [i, j-i+1]
  fill!(u, 0)
  fill!(l, 0)
  # L与U，二者不可或缺
  for i = 1:n
    for j = i:min(i + q, n)
      u[i, j-i+1] = a[i, j]
      for k = max(i - p, j - q, 1):min(i - 1, j - 1)
        u[i, j-i+1] -= l[i-k, k] * u[k, j-k+1]
      end
    end

    # 对称矩阵的福利: L = U' D^-1
    for i2 = i+1:min(i + p, i + q, n)
      l[i2-i, i] = u[i, i2-i+1] / u[i, 1] # u[i, i2]
    end
  end
  l, u
end

function zip_U(U::AbstractMatrix, q)
  U2 = variables(:U, 1:n, 1:q+1)
  fill!(U2, 0)
  for i = 1:n, j = i:min(n, q + i)
    # 0 <= j-i <= q ==> i <= j <= q+i
    U2[i, j-i+1] = U[i, j]
  end
  U2
end
