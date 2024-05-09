export solve_L, solve_U

function solve_U(U::AbstractMatrix{T}, b::AbstractArray) where {T}
  n = length(b)
  x = similar(b)
  x[n] = b[n] / U[n, n]

  @inbounds for i = n-1:-1:1
    x[i] = b[i] - sum(U[i, i+1:n] .* x[i+1:n])
    x[i] /= U[i, i]
  end
  x
end

function solve_L(L::AbstractMatrix{T}, b::AbstractArray) where {T}
  n = length(b)
  x = similar(b)
  x[1] = b[1] / L[1, 1]

  @inbounds for i = 2:n
    x[i] = b[i] - sum(L[i, 1:i-1] .* x[1:i-1])
    x[i] /= L[i, i]
  end
  x
end

################################################################################
## 带状矩阵的版本, kong的数据压缩方式
function solve_U(BU::BandedMat{T}, b::AbstractArray) where {T}
  # [i, j] => [i, j - i + 1]     # 对于U
  (; q) = BU
  U = BU.data
  n = length(b)
  x = similar(b)
  x[n] = b[n] / U[n, 1]

  @inbounds for i = n-1:-1:1
    # x[i] = b[i] - sum(U[i, i+1:n] .* x[i+1:n])
    x[i] = b[i] 
    for k = i+1:min(i+q, n)
      x[i] -= U[i, k-i+1] * x[k]
    end
    x[i] /= U[i, 1]
  end
  x
end

function solve_L(BL::BandedMat{T}, b::AbstractArray) where {T}
  # [i, j] => [i, j - i + p + 1] # 对于L
  (; p) = BL
  L = BL.data
  n = length(b)
  x = similar(b)
  x[1] = b[1]

  @inbounds for i = 2:n
    # x[i] = b[i] - sum(L[i, 1:i-1] .* x[1:i-1])
    x[i] = b[i] 
    for j = max(i-p,1):i-1
      x[i] -= L[i, j-i+p+1] * x[j]
    end
    # x[i] /= L[i, p+1]
  end
  x
end

