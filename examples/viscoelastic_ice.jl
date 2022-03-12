using Oceananigans
using Oceananigans.Units
using GLMakie
using IcyAntics
using IcyAntics: time_step!

#####
##### Domain
#####

L = 100kilometers # square

grid = RectilinearGrid(size = (512, 512),
                       x = (0, L),
                       y = (0, L),
                       halo = (3, 3),
                       topology = (Periodic, Periodic, Flat))

#####
##### Barotropic turbulence ocean simulation
#####

Δx = grid.Δxᶜᵃᵃ
ν₄ = 1/day * Δx^4

ocean_model = NonhydrostaticModel(; grid,
                                  timestepper = :RungeKutta3,
                                  advection = WENO5(),
                                  buoyancy = nothing,
                                  tracers = nothing,
                                  closure = HorizontalScalarBiharmonicDiffusivity(ν=ν₄))

ϵ(x, y, z) = randn()
set!(ocean_model, u=ϵ, v=ϵ)

ocean_simulation = Simulation(ocean_model, Δt=1minute, stop_time=12hours)
run!(ocean_simulation) # generate smooth ocean velocities
ocean_simulation.stop_time = Inf

#####
##### Ice model
#####

ocean_state = (; u=ocean_model.velocities.u, v=ocean_model.velocities.v, ρ=1024.0)
rheology = ViscoElasticRheology(modulus=1.0, viscosity=1e4)
@show ice_Δt = 1e-2 * min(Δx / sqrt(rheology.modulus), Δx^2 / (2 * rheology.viscosity))

ice_model = ContinuumIceModel(; grid, ocean=ocean_state, rheology)
ice_simulation = Simulation(ice_model, Δt=ice_Δt, stop_time=1hour)
run!(ice_simulation)
ice_simulation.stop_time = Inf

@show step_ratio = ceil(Int, ocean_simulation.Δt / ice_simulation.Δt)
ice_simulation.Δt = ocean_simulation.Δt / step_ratio

# A coupled simulation
function coupled_time_step!(ocean_simulation, ice_simulation)
    # Substep ice model
    for substep = 1:step_ratio
        time_step!(ice_simulation)
    end
    time_step!(ocean_simulation)
    return nothing
end

#####
##### Visualization
#####

n = Observable(1) # for visualization

# Ocean speed
uo, vo, wo = ocean_model.velocities
ocean_vorticity = Field(∂x(vo) - ∂y(uo))
ocean_vorticity_n = @lift ($n; interior(compute!(ocean_vorticity))[:, :, 1])

# Ice speed
ui, vi = ice_model.velocities
ice_speed = Field(sqrt(ui^2 + vi^2))
ice_vorticity = Field(∂x(vi) - ∂y(ui))
ice_speed_n = @lift ($n; interior(compute!(ice_speed))[:, :, 1])

# Make figure
x, y, z = nodes((Center, Center, Center), grid)
fig = Figure(resolution=(800, 400))
ax_i = Axis(fig[1, 1], aspect=:equal, xlabel="x (km)", ylabel="y (km)", title="Ice speed")
ax_o = Axis(fig[1, 2], aspect=:equal, xlabel="x (km)", ylabel="y (km)", title="Ocean vorticity")
hm_i = heatmap!(ax_i, 1e-3 * x, 1e-3 * y, ice_speed_n)
hm_o = heatmap!(ax_o, 1e-3 * x, 1e-3 * y, ocean_vorticity_n, colormap=:redblue)
title = @lift ($n; "Coupled ice-ocean simulation after " * prettytime(ice_model.clock.time))
Label(fig[0, :], title)

display(fig)

# Step forward and animate
N = 20
record(fig, "viscoelastic_ice.gif", 1:N, framerate=12) do i
    @info "Plotting frame $i of $N..."
    @time for _ = 1:10
        coupled_time_step!(ocean_simulation, ice_simulation)
    end
    n[] = 1
end

