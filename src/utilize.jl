# funcs = [:size, :-, :+]
# for func in funcs
#   @eval Base.$func(b::BandMat) = Base.$func(b.A)
# end
# Base.getindex(x::BandMat, i, j) = x.A[i, j]
# Base.setindex!(x::BandMat, v, i, j) = x.A[i, j] = v;

function Base.show(io::IO, x::AbstractBandMat{T}) where {T<:Real}
  p = hasfield(typeof(x), :p) ? x.p : 0
  q = hasfield(typeof(x), :q) ? x.q : 0
  zipped = hasfield(typeof(x), :zipped) ? "$(x.zipped)" : "nan"
  printstyled(io, "$(typeof(x)), size=$(x.size), bandwidth=($p, $q), zipped=$zipped\n",
    color=:green, underline=true)
  display(x.data)
  return nothing
end

bandwidth(x::AbstractBandMat) = (x.p, x.q)
bandwidth(x::BandedL) = (x.p, 0)
export bandwidth

## 添加一些简单的函数
ops = [:+, :-, :*, :/]
for fun in ops
  @eval Base.$fun(b::AbstractBandMat, x::Real) = Banded($fun.(b.data, x), b)
  @eval Base.$fun(x::Real, b::AbstractBandMat) = Banded($fun.(x, b.data), b)

  @eval Base.$fun(x::AbstractBandMat, y::AbstractBandMat) = begin
    bandwidth(x) == bandwidth(y) || throw(ArgumentError("bandwidth not match"))
    # if bandwidth(x) == bandwidth(y)
    # band = typeof(x)
    BandedL($fun.(x.data, y.data), x)
  end
end

## lapack
# `A[i, j] => B[j, i - j + q + 1]`

## kong
# A[i, j] => B[i, j - i + p + 1]
# A[j, i] => B[j, i - j + q + 1]
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

# 条带以外的元素填充为0
# force_band!(B::AbstractBandMat) = force_band!(B.data, B.p, B.q)
force_band!(B::AbstractMatrix, p::Int, q::Int) = force_band!(B; p, q)
function force_band!(A::AbstractMatrix{T}; p::Int, q::Int) where {T}
  n, m = size(A)
  for i = 1:n
    for j = 1:min(i - p - 1, m) # 下三角
      A[i, j] = 0
    end
    for j = i+q+1:m # 上三角
      A[i, j] = 0
    end
  end
  return A
end

force_upper!(x::AbstractMatrix) = force_band!(x, 0, size(x, 2) - 1)
force_lower!(x::AbstractMatrix) = force_band!(x, size(x, 1) - 1, 0)

# 强制修改为对称
function force_sym!(A::AbstractMatrix{T}) where {T}
  n, m = size(A)
  @assert n == m "check_sym: matrix is not square"
  @inbounds for i = 1:n
    for j = 1:i-1
      A[i, j] = A[j, i]
    end
  end
  return A
end
