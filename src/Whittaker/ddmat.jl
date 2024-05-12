export ddmat_full, ddmat_band

## 改成带状矩阵的格式
function ddmat_band(x::AbstractVector{T}, d::Integer=2) where {T}
  n = length(x)
  if d == 0
    return BandedMat(ones(T, n, 1), 0, 0; size=(n, n))
  else
    R = diff(ddmat_band(x, d - 1))
    (; q) = R

    for i = 1:n-d
      dx = x[i+d] - x[i]
      for j = 1:q+1
        R.data[i, j] /= dx
      end
    end
    R
  end
end

function ddmat_full(x::AbstractVector{T}, d::Integer=2) where {T}
  n = length(x)
  if d == 0
    return _diagm(ones(T, n))
  else
    R = diff(ddmat_full(x, d - 1))
    m = size(R, 2)

    for i = 1:n-d
      dx = x[i+d] - x[i]
      for j = i:min(i + d, m)
        R[i, j] /= dx
      end
    end
    R
  end
end


function _diagm(x::AbstractVector{T}) where {T}
  n = length(x)
  D = zeros(T, n, n)
  for i = 1:n
    D[i, i] = x[i]
  end
  D
end
