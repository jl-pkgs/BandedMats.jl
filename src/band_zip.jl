band_zip(b::BandMat; kw...) = band_zip(b.data, b.p, b.q; kw...)

function band_zip(A::AbstractMatrix{T}, p::Int, q::Int; type="kong") where {T}
  n, m = size(A)
  # B = zeros(T, p + q + 1, n)
  B = zeros(T, n, p + q + 1)

  if type == "lapack"
    @inbounds for i = 1:n
      for j = max(i - p, 1):min(i + q, m)
        # 0 <= i-j+q <= p+q  ==> i - p <= j <= i + q
        B[j, i-j+q+1] = A[i, j]
        # B[i-j+q+1, j] = A[i, j]
      end
    end
  elseif type == "kong"
    # Whittaker的存储方案
    for i = 1:n
      for k = max(-p, 1 - i):min(q, m - i)
        # 1 <= i + k <= m ==> 1 - i <= k <= m - i
        B[i, k+p+1] = A[i, i+k]
        # j = i+k
        # A[i, j] => B[i, j - i + p + 1]
      end
    end
  end
  return B
end

Base.Matrix(B::BandMat) = B.data

function Base.Matrix(BD::Union{BandedL{T},BandedMat{T}}) where {T}
  p, q = bandwidth(BD)
  (; type) = BD
  n, m = BD.size
  B = BD.data
  A = zeros(T, n, m)

  if type == "lapack"
    @inbounds for i = 1:n
      for j = max(i - p, 1):min(i + q, m)
        # 0 <= i-j+q <= p+q  ==> i - p <= j <= i + q
        # B[j, i-j+q+1] = A[i, j]
        # A[i, j] = B[i-j+q+1, j]
        A[i, j] = B[j, i-j+q+1]
      end
    end
  elseif type == "kong"
    @inbounds for i = 1:n
      for k = max(-p, 1 - i):min(q, m - i)
        A[i, i+k] = B[i, k+p+1]
      end
    end
  end
  return A
end
