# <https://blog.csdn.net/qq_39400324/article/details/123457380>

################################################################################
## 带状矩阵的版本, kong的数据压缩方式
"""
    LDL_solve(BL::BandedMat{T}, d::AbstractArray, b::AbstractArray)

```math
Ax = b
L D L' x = b
L (D L' x) = b
θ = D L' x, θ = L^-1 b
x = L'^-1 D^-1 θ
```

# Examples
```julia
LDL_solve(BL, d, b)
```
"""
function LDL_solve(BL::BandedMat{T}, d::AbstractVector{T}, b::AbstractArray) where {T}
  # [i, j] => [i, j - i + p + 1] # L, A
  # [i, j] => [i, j - i + 1]     # U
  (; p, q) = BL
  L = BL.data
  n = length(b)
  x = similar(b)

  ## L
  x[1] = b[1]
  @inbounds for i = 2:n
    # x[i] = b[i] - sum(L[i, 1:i-1] .* x[1:i-1])
    x[i] = b[i]
    for j = max(i - p, 1):i-1
      x[i] -= L[i, j-i+p+1] * x[j]
    end
  end

  # U的计算
  x[n] = x[n] / d[n]
  @inbounds for i = n-1:-1:1
    x[i] /= d[i]
    for k = i+1:min(i+p, n)
      x[i] -= L[k, i-k+p+1] * x[k] # U[i, k]
      # x[i] -= U[i, k-i+1] * x[k]
    end
  end
  x
end


export LDL_solve, LDL
