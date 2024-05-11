using BenchmarkTools
using BandedMats
using Test

n = 100
A = rand(n, n)
p = q = 4
force_band!(A, p, q)
force_sym!(A)
b = rand(n)

A2 = BandedL(A, p; zipped=false)
A_sym = SymBanded(A2)
@time r_sym = A_sym \ b;
@time r = A \ b;
@test r ≈ r_sym

@btime r_sym = $A_sym \ $b;
# # 60.200 μs (13 allocations: 47.55 KiB)
@btime r = $A \ $b;
# # 35.107 ms (4 allocations: 7.64 MiB)
