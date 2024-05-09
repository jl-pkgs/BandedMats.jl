# include("main_symb.jl")

# funcs = [:size, :-, :+]
# for func in funcs
#   @eval Base.$func(b::BandMat) = Base.$func(b.A)
# end
# Base.getindex(x::BandMat, i, j) = x.A[i, j]
# Base.setindex!(x::BandMat, v, i, j) = x.A[i, j] = v;

function Base.show(io::IO, x::AbstractBandMat{T}) where {T<:Real}
  printstyled(io, "$(typeof(x)): p = $(x.p), q = $(x.q) \n", color=:blue, underline=true)
  display(x.data)
  return nothing
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
