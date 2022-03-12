module IcyAntics

export ContinuumIceModel, ViscoElasticRheology

include("ContinuumIceModels.jl")

using .ContinuumIceModels

end # module
