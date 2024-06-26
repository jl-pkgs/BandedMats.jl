using LinearAlgebra
using BandedMats
using Test



@testset "LU_gauss and LU_full" begin
  function test_lu(size)
    A = rand(size...)
    l, u = lu(A, NoPivot())

    L, U = LU_gauss(A)
    mat_equal(L, l)
    mat_equal(U, u)

    L, U = LU_full(A)
    mat_equal(L, l)
    mat_equal(U, u)
  end
  test_lu((10, 10))
  test_lu((10, 6))
  test_lu((6, 10))
end


@testset "LU_band" begin
  function test_LU_band(A; p, q)
    force_band!(A, p, q)
    # 标准答案
    l, u = lu(A, NoPivot())
    u2 = band_zip(u, 0, q; type="kong")
    l2 = band_zip(l, p, 0; type="kong")[:, 1:end-1]

    # # 我的版本
    B = BandedMat(A, p, q; type="kong", zipped=false)
    # BD = BandedMat(B; type="kong")
    BL, BU = LU_band(B)
    @test BU.data ≈ u2
    @test BL.data ≈ l2

    _l, _u = LU_band_full(A; p, q)
    @test _l ≈ l
    @test _u ≈ u
  end

  n = 6
  test_LU_band(rand(n, n + 2); p=1, q=2)
  test_LU_band(rand(n + 2, n); p=2, q=2)

  test_LU_band(rand(n, n+2); p=2, q=1)
  test_LU_band(rand(n+2, n); p=2, q=1)
  
  test_LU_band(rand(n, n); p=2, q=3)
end
