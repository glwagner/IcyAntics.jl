# IcyAntics.jl

An open-source, community-driven finite-volume sea ice modeling package in the [Oceananigans](https://github.com/CliMA/Oceananigans.jl) ecosystem.

<p align="center">
  <img width=800 src=https://user-images.githubusercontent.com/15271942/158024239-926769a1-eebc-4da3-93ad-9b908c740270.gif>
</p>

_Toy viscoelastic sea ice model. See [`examples/viscoelastic_ice.jl`](https://github.com/glwagner/IcyAntics.jl/blob/main/examples/viscoelastic_ice.jl) or code below._

# Proposal

We propose to develop `IcyAntics.jl` as a sea ice modeling package for both off-line sea ice simulations with prescribed atmosphere-ocean thermomechanical states, and for coupled ice-ocean and ice-atmosphere-ocean simulations at kilometer to planetary scales.

## Mission

`IcyAntics.jl` will enable

* Data-driven sea ice parameterization and model development;
* [Parameter estimation and uncertainty quantification](https://github.com/CliMA/OceanLearning.jl) for new and existing sea ice models;
* High resolution, high performance coupled and uncopuled sea ice simulations with a high-level, user-friendly `Oceananigans`-like interface.

## Core capabilities

Core capabilities will rest on existing, tested sea ice modeling paradigms.
We propose independent efforts to implement the following core capabilities:

* `ContinuumIceModel` with `ElastoViscoPlastic` rheology with either explicit time-stepping or implicit treatment of pseudoelastic dynamics ([Hunke and Dukowicz 1997](https://journals.ametsoc.org/view/journals/phoc/27/9/1520-0485_1997_027_1849_aevpmf_2.0.co_2.xml), [Koldunov et al. 2019](https://agupubs.onlinelibrary.wiley.com/doi/full/10.1029/2018MS001485))
* `DiscreteElementIceModel` Lagrangian ice model ([Chen et al. 2021](https://agupubs.onlinelibrary.wiley.com/doi/full/10.1029/2021MS002513))
* `ContinuumIceModel` with `ViscoPlastic` rheology and implicit time-stepping with a nonlinear implicit solver ([Zhang and Hibler 1997](https://agupubs.onlinelibrary.wiley.com/doi/epdf/10.1029/96JC03744), [Losch et al. 2010](https://www.sciencedirect.com/science/article/pii/S1463500309002418?casa_token=7X6zEGzN43EAAAAA:A1PtPqOSnE-8u9aHyvc2rfffv48yv7sJIbAwyhD1PHb3U_rNFcepGKOMa12wMXqXsI5QDlh4zg))
* `ContinuumIceModel` with `ElastoBrittle` rheology ([Bouillon and Rampal 2015](https://pdf.sciencedirectassets.com/272136/1-s2.0-S1463500315X00060/1-s2.0-S1463500315000694/main.pdf?X-Amz-Security-Token=IQoJb3JpZ2luX2VjELL%2F%2F%2F%2F%2F%2F%2F%2F%2F%2FwEaCXVzLWVhc3QtMSJIMEYCIQCs1kqp%2BBpXiVG60koGCe0Wpc292eMARVIcVDSUXhw7FwIhAJJNoDIOFwDcJD8q%2BCAD1UFXMolOLZ8TibWutXoY7RBmKvoDCCsQBBoMMDU5MDAzNTQ2ODY1IgwBqcjhBxfwG4nLmDoq1wMCAl1spG4A5IMYqlezqs7QsTEB1P0TABxb0yPkKmiwHuVfEYML7I7EfewmwNyqwhWWd0C4fc2nsnpeIB7E%2BKju2ihudPyL70YBSpY50oMqE6QX9Qnt07o8upoOgPXIT%2FND6Qo%2Fmtk4BRcd2uODp6odvSQmPfR5VylC4xxUAyV5W66T0ua7z8Bi0cokn3dSg3ku16J%2F4bnfr2xqJ35nq86iLEl5Q3z%2B94%2FIvzOiyqSByalu6Jaqx3ULHmYOZJwQXTkL7oybTP4Z3o%2BDcl2cPQqDzBBWhYsYW9BLADK3EZr9p9s36JvOxOLcT29zAy3lcK%2FJlQC4WRhFaIbOW6txqYx56XXs9v%2FlJaZtLAag9fpQaXdj0kg6wUHLLqrFNCd83%2BNJvpP7E0afLeMB9fx%2BZ1nR2fjP50sAN2i3U1e4mVE7d8KAinyGholy0uecCMx1IqiHKLoUhvnPV6dHprKkOZCSW0rK51U7KEY4nPQ2PmP15whd%2BwCY4ow5nfFOAXREBVMmhDOpNICGhkVOB8EbUH%2B%2F9AQYCg39e4nTmzq7j%2BxrzquO2rO7JlwvM56IJs8aMQtxRjzknBaG1kQdiT2fNyNDCyCDsfV33RopTHBUta%2BlCtj3uGDbN7kwodaxkQY6pAG0j%2FAzpjaW3Myr97MOPKqQ7a0au2%2BRKG2cUvKO1pYmnZRkq06C8UkWs1QbxF9yhp8zelr3su3cCt13%2B5wLTv1NgKJYYgOI719uRJKKTuFmuHRirTe3HbyP01vLJjUmKGtcuEYd2H6HfOP1WSRWhx9rC3E88wBtFAo0KtCP2g%2B9eFQxziGwXKO7%2FnxAYN2WMnDfUXhJ0yNx68TN5RLAR%2BGnryl9Gw%3D%3D&X-Amz-Algorithm=AWS4-HMAC-SHA256&X-Amz-Date=20220312T111726Z&X-Amz-SignedHeaders=host&X-Amz-Expires=300&X-Amz-Credential=ASIAQ3PHCVTYV7OWEDV7%2F20220312%2Fus-east-1%2Fs3%2Faws4_request&X-Amz-Signature=0ee5b17fa52b239764295cf0f9f6c9ba561fba76df731ce896ab0e75efb0411a&hash=80c54164ef9f7cb37b51fa64a5a89c8c4df85c25a9b275be20b4d8f3d0ab99ce&host=68042c943591013ac2b2430a89b270f6af2c76d8dfd086a07176afe7c76c2c61&pii=S1463500315000694&tid=spdf-3d9100a9-c727-4d95-9c28-4024733711c3&sid=3d8f14755775c04f054a795-b8f297429dffgxrqa&type=client&ua=4c0a030c5001565b5003&rr=6eac1f996c8c314f)) with high-order WENO Eulerian ice advection schemes

## Research directions

Research directions are community-defined.
We _propose_ the following:

* Offline parameter estimation of free parameters in sea ice viscoplastic rheologies, elasto-brittle rheologies, damage models, and therodynamic models;
* Coupled ice-ocean parameter estimation for sea ice models and ocean turbulence parameterizations;
* Extension of discrete element methods to approximate viscoplastic, elasto-brittle, and other sea ice rheology paradigms;
* Numerical methods for sea ice modeling including high-order finite volume formulations and semi-Lagrangian finite volume advection schemes.

## Example code

The following code from [`examples/viscoelastic_ice.jl`](https://github.com/glwagner/IcyAntics.jl/blob/main/examples/viscoelastic_ice.jl) implements a numerical simulation of a "toy" viscoelastic sea ice model with uniform concentration coupled to a barotropic turbulence ocean model:

```julia
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
                                  advection = WENO(),
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

```

## Running example code

1. [Download julia](https://julialang.org/downloads/).

With julia 1.8:

2. git clone https://github.com/glwagner/IcyAntics.jl.git
3. Navigate to `/examples`
4. Type

```
julia --project
```

to open julia in the IcyAntics environment.

5. Type 

```julia
julia> using Pkg; Pkg.instantiate()
```

to install all required packages.

6. Type 

```julia
julia> include("viscoelastic_ice.jl")
```

to run the example.

