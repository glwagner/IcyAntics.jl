module ContinuumIceModels

export ContinuumIceModel, time_step!

# using KernelAbstractions
using Oceananigans.Advection: CenteredSecondOrder
using Oceananigans.BoundaryConditions: fill_halo_regions!
using Oceananigans.Models: AbstractModel
using Oceananigans.AbstractOperations: ∂x, ∂y
using Oceananigans.Fields: XFaceField, YFaceField, CenterField, Field, Center, Face
using Oceananigans.TimeSteppers: Clock, tick!

import Oceananigans.TimeSteppers: time_step!

mutable struct ContinuumIceModel{Grid,
                                 Clo,
                                 Adv,
                                 Rheo,
                                 Drag,
                                 Vel,
                                 Tracers,
                                 Stresses,
                                 Ocean,
                                 Tend} <: AbstractModel{Nothing}
    grid :: Grid
    clock :: Clo
    tracer_advection :: Adv
    rheology :: Rheo
    drag :: Drag
    velocities :: Vel
    tracers :: Tracers # concentration, thickness, snow
    stresses :: Stresses
    ocean :: Ocean # ocean.u, ocean.v, ocean.ρ
    tendencies :: Tend
end

function ContinuumIceModel(; grid, ocean,
                           clock = Clock{eltype(grid)}(0, 0, 1),
                           tracer_advection = CenteredSecondOrder(),
                           rheology = ViscoElasticRheology(1.0, 1.0))

    # Check that grid is 2D with z=0 OR use 3D grid but place fields at the surface.

    velocities = (u = XFaceField(grid), v = YFaceField(grid))

    # One tracer for now
    tracers = (concentration = CenterField(grid), thickness=1)

    stresses = (σ₁₁ = CenterField(grid),
                σ₁₂ = Field{Face, Face, Center}(grid),
                σ₂₂ = CenterField(grid))

    drag = (; Cᴰ_ocean = 1e-3)

    tendencies = (
        u = XFaceField(grid),
        v = YFaceField(grid),
        σ₁₁ = CenterField(grid),
        σ₁₂ = Field{Face, Face, Center}(grid),
        σ₂₂ = CenterField(grid)
    )

    return ContinuumIceModel(grid,
                             clock,
                             tracer_advection,
                             rheology,
                             drag,
                             velocities,
                             tracers,
                             stresses,
                             ocean,
                             tendencies)
end


struct ViscoElasticRheology{T}
    Youngs_modulus :: T
    viscosity :: T
end

Base.@kwdef struct Hibler97Pressure{T}
    p★ :: T = 2.75e5
    c★ :: T = 20.0
end

function calculate_tendencies!(model)
    # Fields
    u, v = model.velocities
    σ₁₁, σ₁₂, σ₂₂ = model.stresses
    Ru, Rv, Rσ₁₁, Rσ₁₂, Rσ₂₂ = model.tendencies

    # Parameters
    E = model.rheology.Youngs_modulus
    μ = model.rheology.viscosity
    Cᴰₒ = model.drag.Cᴰ_ocean
    ρₒ = model.ocean.ρ
    uₒ = model.ocean.u
    vₒ = model.ocean.v

    ϵ₁₁ = ∂x(u)
    ϵ₂₂ = ∂y(v)
    ϵ₁₂ = (∂y(u) + ∂x(v)) / 2

    Rσ₁₁ .= E * ϵ₁₁ - σ₁₁ / 2μ
    Rσ₁₂ .= E * ϵ₁₂ - σ₁₂ / 2μ
    Rσ₂₂ .= E * ϵ₂₂ - σ₂₂ / 2μ

    stress = ∂x(σ₁₁) + ∂y(σ₁₂)
    uᵣ = uₒ - u
    vᵣ = vₒ - v
    drag = Cᴰₒ * ρₒ * sqrt(uᵣ^2 + vᵣ^2) * uᵣ
    Ru .= stress + drag

    return nothing
end

function time_step!(model, Δt)
    calculate_tendencies!(model)

    # Step forward
    fields = merge(model.velocities, model.stresses)

    for (F, dF) in zip(fields, model.tendencies)
        F .+= Δt * dF
    end

    fill_halo_regions!(fields, model.grid.architecture)
    tick!(model.clock, Δt)

    return nothing
end

end # module

