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
