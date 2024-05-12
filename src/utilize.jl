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
