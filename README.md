# IcyAntics.jl

An open-source, community-driven finite-volume sea ice modeling package in the [Oceananigans](https://github.com/CliMA/Oceananigans.jl) ecosystem.

<p align="center">
  <img width=600 src=https://user-images.githubusercontent.com/15271942/158003809-073b3f31-d58a-4883-af3e-755755215a60.gif>
</p>

_Toy viscoelastic sea ice model. See [`examples/viscoelastic_ice.jl`](https://github.com/glwagner/IcyAntics.jl/blob/main/examples/viscoelastic_ice.jl) or code below._

## The proposal

We propose to develop `IcyAntics.jl` as a sea ice modeling package for both off-line sea ice simulations with prescribed atmosphere-ocean thermomechanical states, and for coupled ice-ocean and ice-atmosphere-ocean simulations at kilometer to planetary scales.

## The purpose

`IcyAntics.jl` will enable

* Data-driven sea ice parameterization and model development;
* [Parameter estimation and uncertainty quantification](https://github.com/CliMA/OceanLearning.jl) for new and existing sea ice models;
* High resolution, high performance coupled and uncopuled sea ice simulations with a high-level, user-friendly `Oceananigans`-like interface.

## Core capabilities

Core capabilities will rest on existing, tested sea ice modeling paradigms.
We propose to develop the following core capabilities:

* `ContinuumIceModel` with `ElastoViscoPlastic` rheology and explicit time-stepping for elastic dynamics ([Hunke and Dukowicz 1997](https://journals.ametsoc.org/view/journals/phoc/27/9/1520-0485_1997_027_1849_aevpmf_2.0.co_2.xml))
* `DiscreteElementIceModel` Lagrangian ice model ([Chen et al. 2021](https://agupubs.onlinelibrary.wiley.com/doi/full/10.1029/2021MS002513))
* `ContinuumIceModel` with `ViscoPlastic` rheology and implicit time-stepping with a nonlinear implicit solver ([Zhang and Hibler 1997](https://agupubs.onlinelibrary.wiley.com/doi/epdf/10.1029/96JC03744), [Losch et al. 2010](https://www.sciencedirect.com/science/article/pii/S1463500309002418?casa_token=7X6zEGzN43EAAAAA:A1PtPqOSnE-8u9aHyvc2rfffv48yv7sJIbAwyhD1PHb3U_rNFcepGKOMa12wMXqXsI5QDlh4zg))
* `ContinuumIceModel` with `ElastoBrittle` rheology ([Bouillon and Rampal 2015](https://pdf.sciencedirectassets.com/272136/1-s2.0-S1463500315X00060/1-s2.0-S1463500315000694/main.pdf?X-Amz-Security-Token=IQoJb3JpZ2luX2VjELL%2F%2F%2F%2F%2F%2F%2F%2F%2F%2FwEaCXVzLWVhc3QtMSJIMEYCIQCs1kqp%2BBpXiVG60koGCe0Wpc292eMARVIcVDSUXhw7FwIhAJJNoDIOFwDcJD8q%2BCAD1UFXMolOLZ8TibWutXoY7RBmKvoDCCsQBBoMMDU5MDAzNTQ2ODY1IgwBqcjhBxfwG4nLmDoq1wMCAl1spG4A5IMYqlezqs7QsTEB1P0TABxb0yPkKmiwHuVfEYML7I7EfewmwNyqwhWWd0C4fc2nsnpeIB7E%2BKju2ihudPyL70YBSpY50oMqE6QX9Qnt07o8upoOgPXIT%2FND6Qo%2Fmtk4BRcd2uODp6odvSQmPfR5VylC4xxUAyV5W66T0ua7z8Bi0cokn3dSg3ku16J%2F4bnfr2xqJ35nq86iLEl5Q3z%2B94%2FIvzOiyqSByalu6Jaqx3ULHmYOZJwQXTkL7oybTP4Z3o%2BDcl2cPQqDzBBWhYsYW9BLADK3EZr9p9s36JvOxOLcT29zAy3lcK%2FJlQC4WRhFaIbOW6txqYx56XXs9v%2FlJaZtLAag9fpQaXdj0kg6wUHLLqrFNCd83%2BNJvpP7E0afLeMB9fx%2BZ1nR2fjP50sAN2i3U1e4mVE7d8KAinyGholy0uecCMx1IqiHKLoUhvnPV6dHprKkOZCSW0rK51U7KEY4nPQ2PmP15whd%2BwCY4ow5nfFOAXREBVMmhDOpNICGhkVOB8EbUH%2B%2F9AQYCg39e4nTmzq7j%2BxrzquO2rO7JlwvM56IJs8aMQtxRjzknBaG1kQdiT2fNyNDCyCDsfV33RopTHBUta%2BlCtj3uGDbN7kwodaxkQY6pAG0j%2FAzpjaW3Myr97MOPKqQ7a0au2%2BRKG2cUvKO1pYmnZRkq06C8UkWs1QbxF9yhp8zelr3su3cCt13%2B5wLTv1NgKJYYgOI719uRJKKTuFmuHRirTe3HbyP01vLJjUmKGtcuEYd2H6HfOP1WSRWhx9rC3E88wBtFAo0KtCP2g%2B9eFQxziGwXKO7%2FnxAYN2WMnDfUXhJ0yNx68TN5RLAR%2BGnryl9Gw%3D%3D&X-Amz-Algorithm=AWS4-HMAC-SHA256&X-Amz-Date=20220312T111726Z&X-Amz-SignedHeaders=host&X-Amz-Expires=300&X-Amz-Credential=ASIAQ3PHCVTYV7OWEDV7%2F20220312%2Fus-east-1%2Fs3%2Faws4_request&X-Amz-Signature=0ee5b17fa52b239764295cf0f9f6c9ba561fba76df731ce896ab0e75efb0411a&hash=80c54164ef9f7cb37b51fa64a5a89c8c4df85c25a9b275be20b4d8f3d0ab99ce&host=68042c943591013ac2b2430a89b270f6af2c76d8dfd086a07176afe7c76c2c61&pii=S1463500315000694&tid=spdf-3d9100a9-c727-4d95-9c28-4024733711c3&sid=3d8f14755775c04f054a795-b8f297429dffgxrqa&type=client&ua=4c0a030c5001565b5003&rr=6eac1f996c8c314f)) with high-order WENO Eulerian ice advection schemes

## Research directions

Research directions are community-defined.
We _propose_ the following:

* Offline parameter estimation of free parameters in sea ice viscoplastic rheologies, elasto-brittle rheologies, damage models, and therodynamic models;
* Coupled ice-ocean parameter estimation for sea ice models and ocean turbulence parameterizations;
* Extension of particle-based methods to approximate viscoplastic, elasto-brittle, and other sea ice rheology paradigms;
* Numerical methods for sea ice modeling including high-order finite volume formulations and semi-Lagrangian finite volume advection schemes.

## Example code

The following implements a numerical simulation of a "toy" viscoelastic sea ice model with uniform concentration coupled to an "ocean model":

```julia
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
```
