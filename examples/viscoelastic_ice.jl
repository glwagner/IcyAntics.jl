using Oceananigans
using Oceananigans.Units
using GLMakie
using IcyAntics
using IcyAntics: time_step!

#####
##### Domain: 2D grid with 100km square extent, ~500m resolution
#####

arch = CPU() # change to GPU() to run coupled simulation on GPU.
Nx = Ny = 256
Lx = Ly = 100kilometers

grid = RectilinearGrid(arch,
                       size = (Nx, Ny),
                       halo = (3, 3),
                       x = (0, Lx),
                       y = (0, Ly),
                       topology = (Periodic, Periodic, Flat))

#####
##### Barotropic turbulence ocean simulation
#####

# Ocean mesoscale closure: biharmonic diffusivity
Δh = grid.Δxᶜᵃᵃ
closure = HorizontalScalarBiharmonicDiffusivity(ν = Δh^4 / 4hour)

ocean_model = NonhydrostaticModel(; grid, closure,
                                  timestepper = :RungeKutta3,
                                  advection = WENO5(),
                                  buoyancy = nothing,
                                  tracers = nothing)

# Spin up the ocean model
ϵ(x, y, z) = randn()
set!(ocean_model, u=ϵ, v=ϵ) # random initial velocities
ocean_simulation = Simulation(ocean_model, Δt=2minutes, stop_time=12hours)
run!(ocean_simulation) # generate smooth ocean velocities
ocean_simulation.stop_time = Inf

#####
##### Simple viscoelastic ice model
#####

ocean_state = (; u=ocean_model.velocities.u, v=ocean_model.velocities.v, ρ=1024.0)
rheology = ViscoElasticRheology(modulus=1.0, viscosity=1e4)
@show ice_Δt = 1e-1 * min(Δh / sqrt(rheology.modulus), Δh^2 / (2 * rheology.viscosity))

# Spin up the ice model
ice_model = ContinuumIceModel(; grid, ocean=ocean_state, rheology)
ice_simulation = Simulation(ice_model, Δt=ice_Δt)

#####
##### Construct coupled simulation with substepping
#####

function coupled_time_step!(ocean_simulation, ice_simulation)
    # Floor ice_simulation Δt so `Nsubsteps` ice steps = one ocean step.
    Nsubsteps = ceil(Int, ocean_simulation.Δt / ice_simulation.Δt)
    ice_simulation.Δt = ocean_simulation.Δt / Nsubsteps

    # Substep ice model
    for substep = 1:Nsubsteps
        time_step!(ice_simulation)
    end

    time_step!(ocean_simulation)

    return nothing
end

#####
##### Visualization and coupled run
#####

n = Observable(1) # for visualization

# Ocean vorticity
uo, vo, wo = ocean_model.velocities
ocean_vorticity = Field(∂x(vo) - ∂y(uo))
ocean_vorticity_n = @lift ($n; interior(compute!(ocean_vorticity))[:, :, 1])

# Ice speed
ui, vi = ice_model.velocities
ice_speed = Field(sqrt(ui^2 + vi^2))
ice_speed_n = @lift ($n; interior(compute!(ice_speed))[:, :, 1])

# Make figure
x, y, z = nodes((Center, Center, Center), grid)
fig = Figure(resolution=(800, 400))
ax_i = Axis(fig[1, 1], aspect=1, xlabel="x (km)", ylabel="y (km)", title="Ice speed")
ax_o = Axis(fig[1, 2], aspect=1, xlabel="x (km)", ylabel="y (km)", title="Ocean vorticity")
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

