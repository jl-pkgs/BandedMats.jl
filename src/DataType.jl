abstract type AbstractBandMat{T} end

Base.@kwdef struct BandMat{T} <: AbstractBandMat{T}
  data::AbstractMatrix{T} # how to check value
  p::Int
  q::Int
  size = Base.size(data)
  function BandMat(data::AbstractMatrix{T}, p::Int, q::Int) where {T}
    force_band!(data, p, q) # 地址可能被修改
    new{T}(data, p, q, size(data))
  end
end
BandMat(data, p, q; kw...) = BandMat(; data, p, q, kw...)

Base.@kwdef struct BandedMat{T} <: AbstractBandMat{T}
  data::AbstractMatrix{T} # B
  p::Int
  q::Int
  size = size(data) # original data size
  type = "kong" # "kong", "lapack"
  zipped::Bool = true

  function BandedMat(data::AbstractMatrix{T}, p::Int, q::Int, size, type, zipped) where {T}
    if !zipped
      # force_band!(data, p, q)
      size = Base.size(data)
      data = band_zip(data, p, q; type)
    end
    new{T}(data, p, q, size, type, zipped)
  end
end
BandedMat(data, p, q; kw...) = BandedMat(; data, p, q, kw...)

function BandedMat(b::BandMat; type="kong") 
  B = band_zip(b.data, b.p, b.q; type)
  BandedMat(B, b.p, b.q; type, size=size(b.data))
end

function BandMat(bd::BandedMat{T}) where {T}
  A = band_unzip(bd)
  BandMat(A, bd.p, bd.q)
end


Base.@kwdef struct BandedL{T} <: AbstractBandMat{T}
  data::AbstractMatrix{T} # B
  p::Int
  type = "kong" # "kong", "lapack"
  zipped = true

  function BandedL(data::AbstractMatrix{T}, p::Int, type, zipped) where {T}
    !zipped && (data = band_zip(data, p, 0; type))
    new{T}(data, p, type, zipped)
  end
end

BandedL(data, p; kw...) = BandedL(; data, p, kw...)

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
    @inbounds for i = 1:n
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

# band_unzip(bd::BandedMat) = band_unzip(bd.data, bd.p, bd.q; type=bd.type)
function band_unzip(BD::BandedMat{T}) where {T}
  # function band_unzip(B::BandedMat{T}) where {T}
  (; p, q, type) = BD
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
