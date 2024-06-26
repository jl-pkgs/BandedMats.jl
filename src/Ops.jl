## lapack
# `A[i, j] => B[j, i - j + q + 1]`

## kong
# A[i, j] => B[i, j - i + p + 1]
# A[j, i] => B[j, i - j + q + 1]
function Base.diff(x::AbstractMatrix, d::Integer=1)
  D = x[2:end, :] .- x[1:end-1, :]
  d >= 2 ? diff(D, d - 1) : D
end

# function Base.diff(x::AbstractMatrix, d::Integer=1)
#   # d <= 1 && return x;
#   n, m = size(x)
#   D = zeros(n - 1, m)
#   for j = 1:m
#     for i = 1:n-1
#       D[i, j] = x[i+1, j] - x[i, j]
#     end
#   end
#   d >= 2 ? diff(D, d - 1) : D
# end

function Base.diff(b::BandedMat{T}, d=1) where {T}
  (; p, q) = b
  x = b.data
  n, m = b.size
  _size = (n - 1, m)
  _size2 = (n - 1, p + q + 2) # q+1
  R = BandedMat(zeros(T, _size2),
    p, q + 1; size=_size, zipped=true)
  D = R.data

  @inbounds for i = 2:n
    j_min = max(i - p, 1)
    j_max = min(i + q - 1, m)

    if j_min >= 2
      D[i-1, j_min-i+p+1] = -x[i-1, j_min-i+p+1]
    end

    for j = max(i - p, 1):j_max
      D[i-1, j-i+p+2] = x[i, j-i+p+1] - x[i-1, j-i+p+2]
    end

    if j_max + 1 <= m
      D[i-1, j_max-i+p+3] = x[i, j_max-i+p+2]
    end
  end

  d >= 2 ? diff(R, d - 1) : R
end


function Base.transpose(x::BandedMat{T}) where {T}
  # B[i, j-i+p+1] -> B[j, i - j + q + 1]
  (; p, q) = x
  data = zeros(T, x.size[2], p + q + 1)
  (_n, _m) = x.size
  for i in 1:_n
    for j in max(i - p, 1):min(i + q, _m)
      data[j, i-j+q+1] = x.data[i, j-i+p+1]
    end
  end
  BandedMat(data, q, p; type=x.type)
end

Base.adjoint(x::BandedMat) = transpose(x)

function Base.:\(A::BandedMat{T}, b::AbstractVector{T}) where {T}
  L, U = LU_band(A)
  solve_U(U, solve_L(L, b))
end

function Base.:*(x::BandMat{T1}, y::BandMat{T2}) where {T1,T2}
  p₁, q₁ = x.p, x.q
  p₂, q₂ = y.p, y.q
  n₁, m = x.size
  m, n₂ = y.size

  R = zeros(promote_type(T1, T2), n₁, n₂)
  @inbounds for i = 1:n₁
    for j = 1:n₂
      k_min = max(min(i - p₁, j - q₂), 1)
      k_max = min(max(i + q₁, j + p₂), m)
      for k = k_min:k_max
        # i - p₁ <= k <= i + q₁
        # j - q₂ <= k <= j + p₂  <== -p₁ <= j - k <= q₂  
        R[i, j] += x.data[i, k] * y.data[k, j]
      end
    end
  end
  BandMat(R, p₁ + p₂, q₁ + q₂)
  # return R
end

function Base.:*(x::BandedMat{T1}, y::BandedMat{T2}) where {T1,T2}
  p₁, q₁ = x.p, x.q
  p₂, q₂ = y.p, y.q
  n₁, m = x.size
  m, n₂ = y.size
  # R = zeros(promote_type(T1, T2), n₁, n₂)
  R = zeros(promote_type(T1, T2), n₁, p₁ + p₂ + q₁ + q₂ + 1)
  p = p₁ + p₂
  q = q₁ + q₂

  @inbounds for i = 1:n₁
    for j = 1:n₂
      k_min = max(i - p₁, j - q₂, 1)
      k_max = min(i + q₁, j + p₂, m)
      for k = k_min:k_max
        # R[i, j] += x.data[i, k-i+p₁+1] * y.data[k, j-k+p₂+1]
        R[i, j-i+p+1] += x.data[i, k-i+p₁+1] * y.data[k, j-k+p₂+1]
        # R[i, j-i+p₁+1] += x.data[i, k-i+p₁+1] * y.data[k, j-k+p₂+1]
      end
    end
  end
  BandedMat(R, p, q; type=x.type, zipped=true, size=(n₁, n₂))
end

"""
    U_sq(x::BandedMat{T})

计算时，只使用了上三角的数据。

R = U' * U，同时将结果保存为下三角
"""
function BandedU_sq!(R::BandedL{T}, x::BandedMat{T}) where {T}
  q::Int = x.q
  p::Int = 0
  n::Int, m::Int = x.size
  A = x.data
  r = R.data
  # R = zeros(T, m, m)

  @inbounds for i = 1:m
    for j = i:min(i + q, m)
      k_min = max(i - q, j - q, 1)
      k_max = min(i, j, n)
      for k = k_min:k_max
        r[j, i-j+q+1] += A[k, i-k+p+1] * A[k, j-k+p+1]
      end
    end
  end
  R
end

function BandedU_sq(x::BandedMat{T}) where {T}
  n, m = x.size
  q = x.q
  R = BandedL(zeros(T, m, q+1), q; size=(m, m))
  BandedU_sq!(R, x)
end


export BandedU_sq
