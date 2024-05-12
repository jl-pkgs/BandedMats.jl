# include("GEN_A.jl")

Base.@kwdef mutable struct IntermBand{T}
  n::Int
  p::Int = 3

  # 对称矩阵
  "D2 = D' * D"
  D2::BandedL{T} = BandedU_sq(ddmat_band(1:n, p))

  "A = W + λ * D' * D"
  A::BandedL{T} = BandedL(zeros(T, n, p + 1), p; size=(n, n))

  "L, d = LDL_band(A)"
  L::BandedL{T} = BandedL(zeros(T, n, p), p; size=(n, n))
  d::Vector{T} = zeros(T, n)

  z::Vector{T} = zeros(T, n)
end


"""
A = W + λ * D' * D
"""
function GEN_A!(A::BandedL,
  x::AbstractVector, w::AbstractVector; λ::Real=2.0, p::Int=3)

  D = ddmat_band(x, p)

  BandedU_sq!(A, D) # D2是一个下三角, L
  data = A.data

  @inbounds for i = 1:length(x)
    for j = 1:p+1
      data[i, j] *= λ
    end
    data[i, p+1] += w[i]
  end
  return A
end

# function LDL_band!(interm::IntermBand{T}) where {T}
#   LDL_band!(interm.L, interm.d, interm.A)
# end

# function LDL_solve!(interm::IntermBand{T}, y::AbstractVector{T}) where {T}
#   LDL_solve!(interm.z, interm.L, interm.d, y)
# end

function whit_band(y::AbstractVector{T}, w, x;
  λ::Real=2.0, p::Int=3, interm=nothing) where {T}

  interm === nothing && (interm = IntermBand{T}(; n=length(y), p))
  @unpack L, d, A, z = interm

  GEN_A!(A, x, w; λ, p)
  # L, d = LDL_band(A)
  # LDL_solve(L, d, y .* w)
  LDL_band!(L, d, A)
  LDL_solve!(z, L, d, y .* w)
  return z
end


function whit_band(y::AbstractVector{T}, w;
  λ::Real=2.0, p::Int=3, interm=nothing) where {T}

  interm === nothing && (interm = IntermBand{T}(; n=length(y), p))
  @unpack D2, L, d, A, z = interm

  "A = W + λ * D' * D"
  data = A.data
  @inbounds for i = 1:length(y)
    for j = 1:p+1
      data[i, j] = λ * D2.data[i, j]
    end
    data[i, p+1] += w[i]
  end

  LDL_band!(L, d, A)
  LDL_solve!(z, L, d, y .* w)
  return z
end


export IntermBand
export whit_band

# _D = zeros(n, d + 1)
# D = diff(diagm(x), d)
# D = BandedMat(D, 0, d; zipped=false)

# # A = W + λ * D' * D
# W = BandedL(diagm(w), d; size=(n, n), zipped=false) #|> Matrix
# A = W + BandedL(λ * D' * D)
