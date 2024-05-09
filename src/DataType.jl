abstract type AbstractBandMatrix{T} end

Base.@kwdef struct BandMatrix{T} <: AbstractBandMatrix{T}
  data::AbstractMatrix{T} # how to check value
  p::Int
  q::Int

  function BandMatrix(data::AbstractMatrix{T}, p::Int, q::Int) where {T}
    check_band!(data, p, q) # 地址可能被修改
    new{T}(data, p, q)
  end
end

Base.@kwdef struct BandedMatrix{T} <: AbstractBandMatrix{T}
  data::AbstractMatrix{T} # B
  p::Int
  q::Int
  type = "kong" # "kong", "lapack"
end

# 相互转换
function BandedMatrix(data, p, q; type="kong", zipped=true)
  !zipped && (data = band_zip(data, p, q; type))
  BandedMatrix(; data, p, q, type)
end

function BandedMatrix(b::BandMatrix; type="kong")
  B = band_zip(b.data, b.p, b.q; type)
  BandedMatrix(B, b.p, b.q, type)
end

function BandMatrix(bd::BandedMatrix{T}) where {T}
  A = band_unzip(bd)
  BandMatrix(A, bd.p, bd.q)
end


# 条带以外的元素填充为0
function check_band!(A::AbstractMatrix{T}, p::Int, q::Int) where {T}
  n, m = size(A)
  for i = 1:n
    for j = 1:i-p-1 # 下三角
      A[i, j] = 0
    end
    for j = i+q+1:m # 上三角
      A[i, j] = 0
    end
  end
  A
end

band_zip(b::BandMatrix) = band_zip(b.data, b.p, b.q)
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
    @inbounds for i = 1:n
      for k = max(-p, 1 - i):min(q, m - i)
        # 1 <= i + k <= m ==> 1 - i <= k <= m - i
        B[i, k+p+1] = A[i, i+k]
        # j = i+k
        # A[i, j] => B[i, j - i + p + 1]
      end
    end
  end
  B
end

band_unzip(bd::BandedMatrix) = band_unzip(bd.data, bd.p, bd.q; type=bd.type)
function band_unzip(B::AbstractMatrix{T}, p::Int, q::Int; type="kong") where {T}
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
