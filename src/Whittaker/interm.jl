Base.@kwdef mutable struct interm_whit{T}
  n::Int
  z::Vector{T} = zeros(T, n)

  # c: u1, d: v, e: u2
  d::Vector{T} = zeros(T, n) # v

  c::Vector{T} = zeros(T, n) # u1
  e::Vector{T} = zeros(T, n) # u2
  f::Vector{T} = zeros(T, n) # u2

  s0::Vector{T} = zeros(T, n)
  s1::Vector{T} = zeros(T, n)
  s2::Vector{T} = zeros(T, n)
end
