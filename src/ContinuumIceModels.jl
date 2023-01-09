module ContinuumIceModels

export ContinuumIceModel, time_step!, ViscoElasticRheology

using KernelAbstractions

using Oceananigans.AbstractOperations: ∂x, ∂y
using Oceananigans.Advection: CenteredSecondOrder
using Oceananigans.Architectures: device_event, device
using Oceananigans.BoundaryConditions: fill_halo_regions!
using Oceananigans.Fields: XFaceField, YFaceField, CenterField, Field, Center, Face
using Oceananigans.Models: AbstractModel
using Oceananigans.TimeSteppers: Clock, tick!
using Oceananigans.Utils: launch!

using Oceananigans.Operators

# Simulations interface
import Oceananigans: fields
import Oceananigans.TimeSteppers: time_step!, update_state!
import Oceananigans.Simulations: reset!

mutable struct ContinuumIceModel{Grid,
                                 Tim,
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
    timestepper :: Tim # silly Oceananigans
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

# Oceananigans.Simulations interface
fields(m::ContinuumIceModel) = merge(m.velocities, m.stresses)
update_state!(m::ContinuumIceModel) = fill_halo_regions!(fields(m), m.grid.architecture)
reset!(::Nothing) = nothing # silly Oceananigans # ⟨⟨ eyes ⟩⟩

"""
    ContinuumIceModel(; grid, ocean, rheology, kw...)

Return a continuum model for sea ice on an `ocean` state with `rheology` and `grid`.
"""
function ContinuumIceModel(; grid, ocean,
                           clock = Clock{eltype(grid)}(0, 0, 1),
                           tracer_advection = CenteredSecondOrder(),
                           rheology = ViscoElasticRheology(0.25, 1e3))

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
                             nothing,
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

"""
    ViscoElasticRheology(modulus=0.25, viscosity=1e3)

Return `ViscoElasticRheology` with elastic `modulus`
and dynamic `viscosity` representing viscoelastic ice model in which
the ice stress tensor obeys a prognostic viscoelastic
evolution equation.
"""
Base.@kwdef struct ViscoElasticRheology{T}
    modulus :: T = 0.25 # not old but Young
    viscosity :: T = 1e3
end

# Has free parameters...
Base.@kwdef struct Hibler97Pressure{T}
    p★ :: T = 2.75e5
    c★ :: T = 20.0
    e :: T = 2.0 # yield curve axis ratio
end

#####
##### Time-stepping
#####

@inline ϵ₁₁ᶜᶜᶜ(i, j, k, grid, U) = ∂xᶜᶜᶜ(i, j, k, grid, U.u)
@inline ϵ₁₂ᶠᶠᶜ(i, j, k, grid, U) = (∂yᶠᶠᶜ(i, j, k, grid, U.u) + ∂xᶠᶠᶜ(i, j, k, grid, U.v)) / 2
@inline ϵ₂₂ᶜᶜᶜ(i, j, k, grid, U) = ∂yᶜᶜᶜ(i, j, k, grid, U.v)

""" Calculate the right-hand-side of the free surface displacement (η) equation. """
@kernel function _calculate_tendencies!(R, grid, Uᵢ, Uₒ, σ, E, μ, Cᴰₒ, ρₒ)
    σ₁₁, σ₁₂, σ₂₂ = σ
    i, j = @index(Global, NTuple)
    k = grid.Nz

    @inbounds begin
        # Stress
        # calculate_stress_tendency!(R, i, j, k, grid, rheology, Uᵢ)
        R.σ₁₁[i, j, k] = E * (ϵ₁₁ᶜᶜᶜ(i, j, k, grid, Uᵢ) - σ₁₁[i, j, k] / 2μ)
        R.σ₁₂[i, j, k] = E * (ϵ₁₂ᶠᶠᶜ(i, j, k, grid, Uᵢ) - σ₁₂[i, j, k] / 2μ)
        R.σ₂₂[i, j, k] = E * (ϵ₂₂ᶜᶜᶜ(i, j, k, grid, Uᵢ) - σ₂₂[i, j, k] / 2μ)

        # Momentum
        # calculate_momentum_tendency!(R, i, j, k, grid, rheology, Uᵢ, Uₒ, Uₐ)
        x_stress = ∂xᶠᶜᶜ(i, j, k, grid, σ₁₁) + ∂yᶠᶜᶜ(i, j, k, grid, σ₁₂)
        y_stress = ∂yᶜᶠᶜ(i, j, k, grid, σ₂₂) + ∂xᶜᶠᶜ(i, j, k, grid, σ₁₂) 

        ρᵢ = 910.0 # ice density TODO: make settable
        uᵣ = Uₒ.u[i, j, k] - Uᵢ.u[i, j, k]
        vᵣ = Uₒ.v[i, j, k] - Uᵢ.v[i, j, k]
        x_drag = Cᴰₒ * ρₒ / ρᵢ * sqrt(uᵣ^2 + vᵣ^2) * uᵣ
        y_drag = Cᴰₒ * ρₒ / ρᵢ * sqrt(uᵣ^2 + vᵣ^2) * vᵣ

        R.u[i, j, k] = x_stress + x_drag
        R.v[i, j, k] = y_stress + y_drag

        # Tracers...
    end
end

function calculate_tendencies!(R, grid, Uᵢ, Uₒ, σ, E, μ, Cᴰₒ, ρₒ)
    arch = grid.architecture
    event = launch!(arch, grid, :xy,
                    _calculate_tendencies!,
                    R, grid, Uᵢ, Uₒ, σ, E, μ, Cᴰₒ, ρₒ,
                    dependencies = device_event(arch))

    wait(device(arch), event)
    return nothing
end

function calculate_tendencies!(model)
    # Unpack parameters
    E   = model.rheology.modulus
    μ   = model.rheology.viscosity
    Cᴰₒ = model.drag.Cᴰ_ocean
    ρₒ  = model.ocean.ρ

    calculate_tendencies!(model.tendencies,
                          model.grid,
                          model.velocities,
                          model.ocean,
                          model.stresses,
                          E, μ, Cᴰₒ, ρₒ)

    #=
    # Rate of strain tensor
    ϵ = (∂x(uᵢ), (∂y(uᵢ) + ∂x(vᵢ)) / 2, ∂y(vᵢ))
         
    # Relative velocities
    uᵣ = uₒ - uᵢ
    vᵣ = vₒ - vᵢ

    # σ₁₁, σ₁₂, σ₂₂ = σ
    # ϵ₁₁, ϵ₁₂, ϵ₂₂ = ϵ
    # R.σ₁₁ .= E * (ϵ₁₁ - σ₁₁ / 2μ)
    # R.σ₁₂ .= E * (ϵ₁₂ - σ₁₂ / 2μ)
    # R.σ₂₂ .= E * (ϵ₂₂ - σ₂₂ / 2μ)
    =#

    return nothing
end

@kernel function step_fields!(F, grid, dF, Δt)
    i, j = @index(Global, NTuple)
    k = grid.Nz
    @inbounds begin
        F.u[i, j, k]   += Δt * dF.u[i, j, k]
        F.v[i, j, k]   += Δt * dF.v[i, j, k]
        F.σ₁₁[i, j, k] += Δt * dF.σ₁₁[i, j, k]
        F.σ₁₂[i, j, k] += Δt * dF.σ₁₂[i, j, k]
        F.σ₂₂[i, j, k] += Δt * dF.σ₂₂[i, j, k]
    end
end

function time_step!(model, Δt; callbacks=nothing)
    calculate_tendencies!(model)

    arch = model.grid.architecture
    event = launch!(arch, model.grid, :xy,
                    step_fields!,
                    fields(model), model.grid, model.tendencies, Δt,
                    dependencies = device_event(arch))

    wait(device(arch), event)

    update_state!(model)
    tick!(model.clock, Δt)

    return nothing
end

end # module

