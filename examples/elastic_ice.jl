using Oceananigans
using GLMakie
using IcyAntics
using IcyAntics: time_step!

# Generate an ocean flow
grid = RectilinearGrid(size=(128, 128), x=(0, 2π), y=(0, 2π), halo=(3, 3),
                       topology=(Periodic, Periodic, Flat))

model = NonhydrostaticModel(; grid, timestepper = :RungeKutta3,
                            advection = WENO5(),
                            buoyancy = nothing,
                            tracers = nothing,
                            closure = ScalarDiffusivity(ν=1e-5))

ϵ(x, y, z) = randn()
set!(model, u=ϵ, v=ϵ)
run!(Simulation(model, Δt=0.05, stop_iteration=100))

ocean = (; u=model.velocities.u, v=model.velocities.v, ρ=1024.0)

ice_model = ContinuumIceModel(; grid, ocean)

Δt = 1e-3
time_step!(ice_model, Δt)

fig = Figure()
ax = Axis(fig[1, 1])
n = Observable(1)

u, v = ice_model.velocities
ice_speed = Field(sqrt(u^2 + v^2))

ice_speed_n = @lift begin
    $n
    compute!(ice_speed)
    interior(ice_speed)[:, :, 1]
end

heatmap!(ax, ice_speed_n)

