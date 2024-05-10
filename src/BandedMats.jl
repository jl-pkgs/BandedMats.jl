module BandedMats

# using LinearAlgebra: det
# using SparseArrays

export BandMat, BandedMat, BandedL
export Banded, Band
export band_zip
export force_band!, force_sym!
export force_lower!, force_upper!

include("DataType.jl")
include("utilize.jl")
include("LU.jl")
include("LU_band.jl")
include("LU_solve.jl")
include("LDL_solve.jl")
include("LDL.jl")

include("tools.jl")
include("inv_diag.jl")

Banded = BandedMat
Band = BandMat

end
