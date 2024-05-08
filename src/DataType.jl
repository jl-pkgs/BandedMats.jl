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
  type = "lapack" #, "kong"
end

# 相互转换
BandedMatrix(data, p, q; type="lapack") = BandedMatrix(; data, p, q, type)

function BandedMatrix(b::BandMatrix; type="lapack")
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
