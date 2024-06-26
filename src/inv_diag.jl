"""
  retrieve the diagonal of the inverse of a banded and symmetric matrix B

一种`带状对称矩阵`的快速算法，用于`Whittaker smoother`求解。注意B需要是对称。

根据U矩阵，计算B^-1的对角线元素。

为配合LDL'分解，这里输入变量采用的是L。

```math
B = (U' * D * U)                        # Hutchinson 1985, Eq. 3.1
B^(-1) = B * U^(-1)' + (1 - U) * B^(-1) # Hutchinson 1985, Eq. 3.3
```

# Arguments
- `U`: [
  1 c₁ e₁ f₁ 0
  0 1  c₂ e₂ f₂
  0 0  1  c₃ e₃
  0 0  0  1  c₄
  0 0  0  0  1
]
- `U2`: [
  c₁ e₁ f₁
  c₂ e₂ f₂
  c₃ e₃ 0
  c₄ 0  0
  0  0  0
]

# Examples
```julia
S = inv_diag(L, d)
```
> Dongdong Kong, CUG, 2024-05-07

# References
1. Hutchinson, Michael F., and Frank R. de Hoog. "Smoothing noisy data with spline 
    functions." Numerische Mathematik 47 (1985): 99-106.
"""
function inv_diag(BL::Union{BandedMat{T},BandedL{T}}, d::AbstractVector{T}) where {T<:Real}
  L::Matrix{T} = BL.data # [n, m]
  m = BL.p
  n = size(L, 1)

  B = zeros(T, n, m + 1)
  B[n, 1] = 1 / d[n]

  @inbounds for i = n-1:-1:1
    B[i, 1] = 1 / d[i]
    for l = 1:min(m, n - i)
      B[i, 1+l] = 0
      for k = 1:min(m, n - i)
        # k <= l, u[i, i+k] * b[i+k, i+l] => u[i, k+1] * b[i+k, l-k+1]
        # k >  l, u[i, i+k] * b[i+l, i+k] => u[i, k+1] * b[i+l, k-l+1]
        _i, _j = k <= l ? (i + k, l - k + 1) : (i + l, k - l + 1)

        B[i, 1+l] -= L[i+k, m+1-k] * B[_i, _j]
        # println("(i=$i, l=$l, k=$k) = $(L[i+k, m+1-k])")
        # B[i, 1+l] -= U[i, k] * B[_i, _j]
      end

      B[i, 1] -= L[i+l, m+1-l] * B[i, 1+l]
      # B[i, 1] -= U[i, l] * B[i, 1+l]
    end
  end
  return B[:, 1]
end

"""
A = U' D U (A = L D L')
U = L' # ! note
"""
function cal_diag(U2::AbstractMatrix{T}, d::AbstractVector{T}; m=3) where {T<:Real}
  n = length(d)
  B = zeros(T, n, m + 1)
  B[n, 1] = 1 / d[n]

  for i = n-1:-1:1
    B[i, 1] = 1 / d[i]
    for l = 1:min(m, n - i)
      B[i, 1+l] = 0
      for k = 1:min(n - i, m)
        _i, _j = k <= l ? (i + k, l - k + 1) : (i + l, k - l + 1)
        B[i, 1+l] -= U2[i, k] * B[_i, _j]
      end
      B[i, 1] -= U2[i, l] * B[i, 1+l]
    end
  end
  B[:, 1]
end

"""
A = L D L'
U = D L'
"""
function inv_diag(U2::AbstractMatrix, d::AbstractVector{T}; m::Int) where {T<:Real}
  n = size(U2, 1)
  B = zeros(T, n, m + 1)
  B[n, 1] = 1 / d[n]

  @inbounds for i = n-1:-1:1
    B[i, 1] = 1.0
    for l = 1:min(m, n - i)
      B[i, 1+l] = 0
      for k = 1:min(m, n - i)
        # k <= l, u[i, i+k] * b[i+k, i+l] => u[i, k+1] * b[i+k, l-k+1]
        # k >  l, u[i, i+k] * b[i+l, i+k] => u[i, k+1] * b[i+l, k-l+1]
        _i, _j = k <= l ? (i + k, l - k + 1) : (i + l, k - l + 1)

        # B[i, 1+l] -= L[i+k, m+1-k] * B[_i, _j]
        B[i, 1+l] -= U2[i, k] * B[_i, _j]
      end
      B[i, 1+l] /= d[i]

      # B[i, 1] -= L[i+l, m+1-l] * B[i, 1+l]
      B[i, 1] -= U2[i, l] * B[i, 1+l]
    end
    B[i, 1] /= d[i]
  end
  return B[:, 1]
end

export inv_diag, cal_diag
