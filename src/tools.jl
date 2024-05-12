function Base.diff(x::AbstractMatrix, d::Integer=1)
  D = x[2:end, :] .- x[1:end-1, :]
  d >= 2 ? diff(D, d - 1) : D
end

# function Base.diff(x::SparseMatrixCSC, d::Integer=1)
#   D = x[2:end, :] .- x[1:end-1, :]
#   d >= 2 ? diff(D, d - 1) : D
# end

# speye(n) = SparseArrays.sparse(I, n, n)

# function ddmat(x::AbstractVector, d::Integer=2)
#   n = length(x)
#   if d == 0
#     return speye(n)
#   else
#     return diff(ddmat(x, d - 1))
#   end
# end



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
