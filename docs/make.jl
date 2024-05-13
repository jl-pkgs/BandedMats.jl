using Documenter, BandedMats
CI = get(ENV, "CI", nothing) == "true"

# Logging.disable_logging(Logging.Warn)

# Make the docs, without running the tests again
# We need to explicitly add all the extensions here
makedocs(
  modules=[
    BandedMats
  ],
  format=Documenter.HTML(
    prettyurls=CI,
  ),
  pages=[
    "Home" => "index.md",
  ],
  sitename="BandedMats.jl",
  warnonly=true,
  clean=false,
)

# Enable logging to console again
# Logging.disable_logging(Logging.BelowMinLevel)
deploydocs(
  repo="github.com/jl-pkgs/BandedMats.jl.git",
)
