# funcs = [:size, :-, :+]
# for func in funcs
#   @eval Base.$func(b::BandMat) = Base.$func(b.A)
# end
# Base.getindex(x::BandMat, i, j) = x.A[i, j]
# Base.setindex!(x::BandMat, v, i, j) = x.A[i, j] = v;

function Base.show(io::IO, x::AbstractBandMat{T}) where {T<:Real}
  p = hasfield(typeof(x), :p) ? x.p : 0
  q = hasfield(typeof(x), :q) ? x.q : 0
  printstyled(io, "$(typeof(x)): p = $p, q = $q \n", color=:blue, underline=true)
  display(x.data)
  return nothing
end

## lapack
# `A[i, j] => B[j, i - j + q + 1]`

## kong
# A[i, j] => B[i, j - i + p + 1]
# A[j, i] => B[j, i - j + q + 1]
function Base.transpose(x::BandedMat{T}) where {T}
  # B[i, j-i+p+1] -> B[j, i - j + q + 1]
  (; p, q) = x
  n, m = size(x.data)
  data = zeros(T, n, m)

  for i = 1:n
    for j = max(1, p - i + 2):min(m, p - i + 2 + m)
      # _i = i
      # _j = j + i - p - 1
      # i2 = _j
      # j2 = _i - _j + q + 1 = m - j + 1
      data[j+i-p-1, m-j+1] = x.data[i, j]
    end
  end
  BandedMat(data, q, p; type=x.type)
end

Base.adjoint(x::BandedMat) = transpose(x)


function Base.:*(x::BandMat{T1}, y::BandMat{T2}) where {T1,T2}
  p₁, q₁ = x.p, x.q
  p₂, q₂ = y.p, y.q
  n₁, m = size(x.data)
  m, n₂ = size(y.data)

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
  return R
end

function Base.:*(x::BandedMat{T1}, y::BandedMat{T2}) where {T1,T2}
  p₁, q₁ = x.p, x.q
  p₂, q₂ = y.p, y.q
  n₁, m = size(x.data)
  m, n₂ = size(y.data)

  R = zeros(promote_type(T1, T2), n₁, n₂)
  @inbounds for i = 1:n₁
    for j = 1:n₂
      k_min = max(min(i - p₁, j - q₂), 1)
      k_max = min(max(i + q₁, j + p₂), m)
      for k = k_min:k_max
        # i - p₁ <= k <= i + q₁
        # j - q₂ <= k <= j + p₂  <== -p₂ <= j - k <= q₂  
        R[i, j] += x.data[i, k-i+p₁+1] * y.data[k, j-k+p₂+1]
      end
    end
  end
  return R
end

# 条带以外的元素填充为0
function force_band!(A::AbstractMatrix{T}, p::Int, q::Int) where {T}
  n, m = size(A)
  @inbounds for i = 1:n
    for j = 1:i-p-1 # 下三角
      A[i, j] = 0
    end
    for j = i+q+1:m # 上三角
      A[i, j] = 0
    end
  end
  return A
end

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

# """
# 代数余子式(algebraic complement)
# """
# function complement(A::AbstractArray, i, j; verbose=false)
#   i, j = j, i
#   m, n = size(A)
#   _A = A[setdiff(1:m, i), setdiff(1:n, j)]
#   verbose && display(_A)
#   (-1)^(i + j) * det(_A)
# end

# function complement(A::AbstractArray{T}) where {T}
#   R = zeros(T, size(A))
#   m, n = size(A)
#   for i = 1:m, j = 1:n
#     R[i, j] = complement(A, i, j) # 代数余子式需要进行一次转置，才能得到A* = C'
#   end
#   R
# end

# export complement
