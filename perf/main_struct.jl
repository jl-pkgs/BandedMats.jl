using Mixers

abstract type AbstractBandMat{T} end

@premix Base.@kwdef struct BAND{T}
  data::AbstractMatrix{T}
  # size::Tuple{Int64,Int64}
  p::Int
end

@premix Base.@kwdef mutable struct BANDED{T}
  data::AbstractMatrix{T}
  size::Tuple{Int64,Int64}
  bandwidth::Tuple{Int64,Int64}
  zipped = true
  storage = "kong"
  # p::Int
  # 添加一个zipped function
end

# 对称矩阵，只保留了L的数据
@BAND struct Band{} end
@BANDED struct Banded{} end
@BANDED struct BandedL{} end
@BANDED struct BandedU{} end
@BANDED struct SymBanded{} end

function BandedL(data::AbstractMatrix{T}, size, bandwidth, zipped=true, storage="kong") where {T}
  p, q = bandwidth
  q == 0 || throw(ArgumentError("BandedL only supports lower triangular matrices"))
  if !zipped
    size = Base.size(data)
    data = band_zip(data, p, q; type)
    zipped = true
  end
  new{T}(data, size, bandwidth, zipped, storage)
end

function BandedU(data::AbstractMatrix{T}, size, bandwidth, zipped=true, storage="kong") where {T}
  p, q = bandwidth
  p == 0 || throw(ArgumentError("BandedL only supports lower triangular matrices"))
  if !zipped
    size = Base.size(data)
    data = band_zip(data, p, q; type)
    zipped = true
  end
  new{T}(data, size, bandwidth, zipped, storage)
end

function SymBanded(data::AbstractMatrix{T}, size, bandwidth, zipped=true, storage="kong") where {T}
  p, q = bandwidth
  p == q || throw(ArgumentError("BandedL only supports lower triangular matrices"))
  if !zipped
    size = Base.size(data)
    # 仅保留下三角的数据
    data = band_zip(data, p, 0; type)
    # force_sym!(data)
    zipped = true
  end
  new{T}(data, size, (p, 0), zipped, storage)
end


# BandedL(rand(4, 4), (4, 4), (2, 0), true, "kong")
BandedL(rand(4, 4), (4, 4), (2, 0))
BandedU(rand(4, 4), (4, 4), (2, 0))
SymBanded(rand(4, 4), (4, 4), (2, 2))
