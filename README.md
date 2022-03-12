# IcyAntics.jl

Under development: a sea ice modeling framework in the [`Oceananigans.jl`](https://github.com/CliMA/Oceananigans.jl) and [`OceanLearning.jl`](https://github.com/CliMA/OceanLearning.jl) ecosystem.
`IcyAntics.jl` is designed to enable parameter estimation and uncertainty quantification for sea ice models, data-driven sea ice parameterization development, and sea ice modeling from kilometer to planetary scales.

Models we might implement:

* Static viscoplastic continuum models ([Zhang and Hibler 1997](https://agupubs.onlinelibrary.wiley.com/doi/epdf/10.1029/96JC03744), [Losch et al. 2010](https://www.sciencedirect.com/science/article/pii/S1463500309002418?casa_token=7X6zEGzN43EAAAAA:A1PtPqOSnE-8u9aHyvc2rfffv48yv7sJIbAwyhD1PHb3U_rNFcepGKOMa12wMXqXsI5QDlh4zg))
* Dynamic elastic-viscous-plastic continuum models ([Hunke and Dukowicz 1997](https://journals.ametsoc.org/view/journals/phoc/27/9/1520-0485_1997_027_1849_aevpmf_2.0.co_2.xml))
* Discrete floe models ([Chen et al. 2021](https://agupubs.onlinelibrary.wiley.com/doi/full/10.1029/2021MS002513))

Sample code (see `examples/elastic_ice.jl`):

```julia
using Oceananigans
using GLMakie
using IcyAntics
using IcyAntics: time_step!

# A simple ocean simulation
grid = RectilinearGrid(size=(128, 128), x=(0, 2π), y=(0, 2π), halo=(3, 3),
                       topology=(Periodic, Periodic, Flat))

ocean_model = NonhydrostaticModel(; grid, timestepper = :RungeKutta3,
                                  advection = WENO5(),
                                  buoyancy = nothing,
                                  tracers = nothing,
                                  closure = ScalarDiffusivity(ν=1e-5))

ϵ(x, y, z) = randn()
set!(ocean_model, u=ϵ, v=ϵ)
ocean_simulation = Simulation(ocean_model, Δt=0.01, stop_iteration=500)
run!(ocean_simulation)
ocean_simulation.stop_iteration = Inf

# A continuum ice model
ocean = (; u=ocean_model.velocities.u, v=ocean_model.velocities.v, ρ=1024.0)
ice_model = ContinuumIceModel(; grid, ocean)

# A coupled model
function coupled_steps!(ocean_simulation, ice_model, N=10)
    for ocean_step = 1:N
        for substep = 1:10
            time_step!(ice_model, 1e-3)
        end
        time_step!(ocean_simulation)
    end
end

# Visualize and step forward:
fig = Figure(resolution=(1200, 600))
ax_i = Axis(fig[1, 1], xlabel="x", ylabel="y", title="Ice speed")
ax_o = Axis(fig[1, 2], xlabel="x", ylabel="y", title="Ocean speed")
n = Observable(1) # for visualization

# Ocean speed visualization
uo, vo, wo = ocean_model.velocities
ocean_speed = Field(sqrt(uo^2 + vo^2))

ocean_speed_n = @lift begin
    $n
    compute!(ocean_speed)
    interior(ocean_speed)[:, :, 1]
end

# Ice speed visualization
ui, vi = ice_model.velocities
ice_speed = Field(sqrt(ui^2 + vi^2))

ice_speed_n = @lift begin
    $n
    compute!(ice_speed)
    interior(ice_speed)[:, :, 1]
end

heatmap!(ax_i, ice_speed_n)
heatmap!(ax_o, ocean_speed_n)

display(fig)

N = 100
record(fig, "elastic_ice.gif", 1:N, framerate=12) do i
    @info "Plotting frame $i of $N..."
    coupled_steps!(ocean_simulation, ice_model)
    n[] = 1
end
```
