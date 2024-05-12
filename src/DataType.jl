abstract type AbstractBandMat{T} end

Base.@kwdef struct BandMat{T} <: AbstractBandMat{T}
  data::Matrix{T} # how to check value
  p::Int
  q::Int
  size::Tuple{Int64,Int64} = Base.size(data)

  function BandMat(data::AbstractMatrix{T}, p::Int, q::Int) where {T}
    force_band!(data, p, q) # 地址可能被修改
    new{T}(data, p, q, size(data))
  end
end

Base.@kwdef struct BandedMat{T} <: AbstractBandMat{T}
  data::Matrix{T} # B
  p::Int
  q::Int
  size::Tuple{Int64,Int64} = size(data) # original data size
  type::String = "kong" # "kong", "lapack"
  zipped::Bool = true

  function BandedMat(data::AbstractMatrix{T}, p::Int, q::Int, size, type, zipped) where {T}
    if !zipped
      # force_band!(data, p, q)
      size = Base.size(data)
      data = band_zip(data, p, q; type)
      zipped = true
    end
    new{T}(data, p, q, size, type, zipped)
  end
end

# abstract type BandedL2{T} <: BandedMat{T} end
Base.@kwdef struct BandedL{T} <: AbstractBandMat{T}
  data::Matrix{T} # B
  p::Int
  size::Tuple{Int64,Int64} = Base.size(data)
  type::String = "kong" # "kong", "lapack"
  zipped::Bool = true

  function BandedL(data::AbstractMatrix{T}, p::Int, size, type, zipped) where {T}
    if !zipped
      size = Base.size(data)
      data = band_zip(data, p, 0; type)
      zipped = true
    end
    new{T}(data, p, size, type, zipped)
  end
end

Base.@kwdef struct SymBandedMat{T} <: AbstractBandMat{T}
  data::Matrix{T} # how to check value
  p::Int
  size::Tuple{Int64,Int64} = Base.size(data)
  type::String = "kong" # "kong", "lapack"
  zipped::Bool = true

  function SymBandedMat(data::AbstractMatrix{T}, p::Int, size, type, zipped) where {T}
    if !zipped
      size = Base.size(data)
      data = band_zip(data, p, 0; type)
      zipped = true
    end
    new{T}(data, p, size, type, zipped)
  end
end




BandMat(data, p, q; kw...) = BandMat(; data, p, q, kw...)

function BandMat(bd::BandedMat{T}) where {T}
  BandMat(Matrix(bd), bd.p, bd.q)
end

BandedMat(data, p, q; kw...) = BandedMat(; data, p, q, kw...)

BandedMat(data, b::BandedMat) = BandedMat(data, b.p, b.q, b.size, b.type, b.zipped)

function BandedMat(b::BandMat; type="kong")
  B = band_zip(b.data, b.p, b.q; type)
  BandedMat(B, b.p, b.q; type, size=size(b.data))
end


BandedL(data::AbstractMatrix, p::Int; kw...) = BandedL(; data, p, kw...)

function BandedL(data::AbstractMatrix, b::BandedL)
  (; p, type, zipped, size) = b
  BandedL(; data, p, type, zipped, size)
end

function BandedL(b::BandedMat)
  (; p, type, zipped, size) = b
  BandedL(b.data[:, 1:1+p], p; type, zipped, size)
end



SymBandedMat(data::AbstractMatrix, p::Int; kw...) = SymBandedMat(; data, p, kw...)

function SymBandedMat(x::BandedL)
  (; p, size, type, zipped) = x
  SymBandedMat(x.data, p; size, type, zipped)
end

SymBanded = SymBandedMat
