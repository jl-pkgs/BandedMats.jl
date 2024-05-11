# 生成D矩阵的构造函数
function coef_diff(d::Int)
  x = zeros(Int, d + 1)
  for i = 1:d+1
    x[i] = (-1)^(i + d + 1) * binomial(d, i - 1)
  end
  return x
end

# 只求下三角的部分
function coef_mult(i::Int, j::Int; n::Int, p::Int, coef::AbstractVector)
  m = p + 1
  k = abs(j - i) # 隔几位
  I = k+1:m
  J = 1:m-k

  k > p && return 0

  if i <= p && j <= p
    num = min(i, j)
    I = I[1:num]
    J = J[1:num]
  elseif i >= n - p + 1 && j >= n - p + 1
    num = n - i + 1
    I = I[end-num+1:end]
    J = J[end-num+1:end]
  end
  
  ans = 0
  for (_i, _j) = zip(I, J)
    ans += coef[_i] * coef[_j]
  end
  return ans
end

## 一步到位A矩阵的构造函数
function GEN_A(x::AbstractVector{T1}, w::AbstractVector{T2}; λ=2.0, d=3) where{T1, T2}
  A = BandedL(zeros(promote_type(T1, T2), n, d + 1), d; size=(n, n))
  GEN_A!(A, x, w; λ, d)
  # BandedL(data, d; size=(n, n))
end

function GEN_A!(A::BandedL{T}, x::AbstractVector{T}, w::AbstractVector{T}; 
  λ=2.0, p::Int=3) where {T}

  n = length(x)
  # data = zeros(n, p + 1)
  data = A.data
  coef::Vector{T} = coef_diff(p)
  λ = T(λ)
  
  @inbounds for i = 1:n
    for j = max(i - p, 1):i
      c = coef_mult(i, j; n, p, coef)
      data[i, j-i+p+1] = x[i] * x[j] * c * λ
    end
    data[i, p+1] += w[i]
  end
  return A
end

Base.@kwdef mutable struct IntermBand{T}
  n::Int
  p::Int = 3
  A::BandedL = BandedL(zeros(T, n, p + 1), p; size=(n, n))
  L::BandedL = BandedL(zeros(T, n, p + 1), p; size=(n, n))
  d::Vector{T} = zeros(T, n)
  z::Vector{T} = zeros(T, n)
end

# import BandedMats: LDL_band!, LDL_solve!

function LDL_band!(interm::IntermBand{T}) where {T}
  LDL_band!(interm.L, interm.d, interm.A)
end

function LDL_solve!(interm::IntermBand{T}, y::AbstractVector{T}) where {T}
  LDL_solve!(interm.z, interm.L, interm.d, y)
end

function whit_band(y::AbstractVector{T}, w, x; λ, p, interm=nothing) where {T}
  interm === nothing && (interm = IntermBand{T}(; n=length(y), p))

  GEN_A!(interm.A, x, w; λ, p)
  LDL_band!(interm)
  LDL_solve!(interm, y.*w)
  interm.z
end

export IntermBand
export whit_band

## A矩阵要保留空间
# _D = zeros(n, d + 1)
# D = diff(diagm(x), d)
# D = BandedMat(D, 0, d; zipped=false)

# # A = W + λ * D' * D
# W = BandedL(diagm(w), d; size=(n, n), zipped=false) #|> Matrix
# A = W + BandedL(λ * D' * D)
