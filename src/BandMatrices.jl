module BandMatrices

# using LinearAlgebra: det
# using SparseArrays

export BandMatrix, BandedMatrix
export band_zip, band_unzip, check_band!

include("DataType.jl")
include("band_zip.jl")
include("utilize.jl")
include("LU.jl")
include("LU_solve.jl")

end
