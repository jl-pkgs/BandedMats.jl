band_zip(b::BandMatrix) = band_zip(b.data, b.p, b.q)
function band_zip(A::AbstractMatrix{T}, p::Int, q::Int; type="lapack") where {T}
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
    @inbounds for i = 1:n
      for k = max(-p, 1 - i):min(q, m - i)
        # 1 <= i + k <= m ==> 1 - i <= k <= m - i
        B[i, k+p+1] = A[i, i+k]
        # (k + p + 1) - (i + k) = p + 1 - i
        # A[i, j] => B[i, p + 1 - i]
      end
    end
  end
  B
end

band_unzip(bd::BandedMatrix) = band_unzip(bd.data, bd.p, bd.q; type=bd.type)
function band_unzip(B::AbstractMatrix{T}, p::Int, q::Int; type="lapack") where {T}
  # function band_unzip(B::BandedMatrix{T}) where {T}
  # (; p, q, type) = B
  n = size(B, 1)
  A = zeros(T, n, n)

  if type == "lapack"
    @inbounds for i = 1:n
      for j = max(i - p, 1):min(i + q, n)
        # 0 <= i-j+q <= p+q  ==> i - p <= j <= i + q
        # B[j, i-j+q+1] = A[i, j]
        # A[i, j] = B[i-j+q+1, j]
        A[i, j] = B[j, i-j+q+1]
      end
    end
  elseif type == "kong"
    @inbounds for i = 1:n
      for k = max(-p, 1 - i):min(q, n - i)
        A[i, i+k] = B[i, k+p+1]
      end
    end
  end
  A
end
