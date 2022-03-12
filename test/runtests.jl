using Test

# eg
#
# $ julia --project
#
# Then
#
# julia> using Pkg; Pkg.test()

icy_antics() = "Is a sea ice modeling framework"
@test icy_antics() == "Is a sea ice modeling framework"
