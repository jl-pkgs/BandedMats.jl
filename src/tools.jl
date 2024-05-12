function copy_L2U!(A)
  for i = axes(A, 1)
    for j = i+1:size(A, 2)
      A[i, j] = A[j, i]
    end
  end
  A
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
