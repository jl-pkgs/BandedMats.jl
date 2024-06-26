using Test
using LinearAlgebra
using SparseArrays
using BenchmarkTools
using BandedMats

# ddmat(x::AbstractVector, d::Integer=2) = diff(diagm(x), d)
ddmat(x::AbstractVector, d::Integer=2) = diff(spdiagm(x), d)

function WHIT(y::AbstractVector, w::AbstractVector; kw...)
  WHIT(y, w, 1:length(y); kw...)
end

function WHIT(y::AbstractVector, w::AbstractVector, x::AbstractVector;
  λ=2.0, p=2)
  n = length(y)
  D = ddmat(x, p)

  W = spdiagm(w)
  A = W + λ * D' * D

  # L = cholesky(A).L # Matrix
  L = cholesky(A, perm=1:n).L # sparse
  z = L' \ (L \ (w .* y))
  z
end


n = 100
y = rand(n)
w = rand(n)
x = rand(n)
λ, p = 2.0, 3

## 测试运行速度
interm = IntermBand{Float64}(; n=length(y), p=3)
# @profview 
# @time

# 时间都浪费在生成D矩阵了
# @time 
@time for i = 1:100_000
  # z1 = whit_band(y, w; λ=2.0, p=3, interm)
  z2 = whit3(y, w; λ=2.0, include_cve=false)
  # z2 = WHIT(y, w, x; λ=2.0, p=3)
end

z1 = whit_band(y, w; λ=2.0, p=3, interm)
z2 = whit3(y, w; λ=2.0, include_cve=false)[1]
z1 - z2

# @btime z1 = whit_band($y, $w, $x; λ=2.0, p=3, interm);
@btime z1 = whit_band($y, $w; λ=2.0, p=3, interm=$interm);
@btime z2 = whit3($y, $w; λ=2.0, include_cve=false);
