export ddmat_full, ddmat_band


function ddmat_band(x::AbstractVector{T}, d::Integer=2) where {T}
  n = length(x)
  B = BandedMat(zeros(n, d + 1), 0, d; size=(n, n))
  R = zeros(n, d + 1)

  fill!(R, 0)
  R[1:n, 1] .= 1

  for k = 1:d
    for i = 1:n-k
      dx = (x[i+k] - x[i]) / k
      R[i, 1] = -R[i, 1] / dx

      for j = i+1:min(i + k, n)
        # R[i, j] = (R[i+1, j] - R[i, j]) / dx
        R[i, j-i+1] = (R[i+1, j-i] - R[i, j-i+1]) / dx
      end
    end
  end

  BandedMat(R[1:n-d, :], 0, d; size=(n - d, n))
  # B
end


function ddmat_full(x::AbstractVector{T}, d::Integer=1) where {T}
  n = length(x)
  R = zeros(T, n, n)
  ddmat_full!(R, x, d)
end

function ddmat_full!(R::AbstractMatrix{T}, x::AbstractVector{T}, d::Integer=1) where {T}
  n = length(x)
  fill!(R, 0)
  @inbounds for i = 1:n
    R[i, i] = 1
  end
  # R = _diagm(ones(T, n))

  m = size(R, 2)
  for k = 1:d
    for i = 1:n-k
      dx = (x[i+k] - x[i]) / k
      for j = i:min(i + k, m)
        R[i, j] = (R[i+1, j] - R[i, j]) / dx
      end
    end
  end
  return @view R[1:n-d, :]
end

# function _diagm(x::AbstractVector{T}) where {T}
#   n = length(x)
#   D = zeros(T, n, n)
#   for i = 1:n
#     D[i, i] = x[i]
#   end
#   D
# end

# function ddmat_full(x::AbstractVector{T}, d::Integer=2) where {T}
#   n = length(x)
#   if d == 0
#     return _diagm(ones(T, n))
#   else
#     R = diff(ddmat_full(x, d - 1))
#     m = size(R, 2)
#     for i = 1:n-d
#       dx = (x[i+d] - x[i]) / d
#       for j = i:min(i + d, m)
#         R[i, j] /= dx
#       end
#     end
#     R
#   end
# end

## 改成带状矩阵的格式
# function ddmat_band(x::AbstractVector{T}, d::Integer=2) where {T}
#   n = length(x)
#   if d == 0
#     return BandedMat(ones(T, n, 1), 0, 0; size=(n, n))
#   else
#     R = diff(ddmat_band(x, d - 1))
#     q::Int = R.q

#     @inbounds for i = 1:n-d
#       dx = (x[i+d] - x[i]) / d
#       for j = 1:q+1
#         R.data[i, j] /= dx
#       end
#     end
#     R
#   end
# end

export ddmat_band2
