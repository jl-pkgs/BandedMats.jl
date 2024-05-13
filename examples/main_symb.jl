using Symbolics
import Symbolics: scalarize, variables
using BandedMats
using LinearAlgebra

@variables λ
Vec(name, n) = variables(name, 1:n)
Mat(name, n) = variables(name, 1:n, 1:n)
Mat_band(name, n; p=2, q=1) = force_band!(Mat(name, n), p, q)

a = Mat(:a, 5)
p, q = 2, 1

x = Vec(:a, 5)
force_lower!(Mat(:a, 5))
force_upper!(Mat(:a, 5))
diagm(x)
