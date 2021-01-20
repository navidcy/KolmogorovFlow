# KolmogorovFlow

Study of the stability of Kolmogorov flow. The code is written in [Julia](https://julialang.org) and utilizes the Julia package [GeophysicalFlows.jl](http://github.com/FourierFlows/GeophysicalFlows.jl).

## Installation

First you need to [install Julia](https://julialang.org/downloads/). Then clone the repository, e.g.,

```
git clone https://github.com/navidcy/KolmogorovFlow.git
```

Enter the directory you've cloned the repository, e.g., 

```
cd KolmogorovFlow
```

Then, after you edit `setup_and_run_simulation.jl` file with your parameter values, and while still inside the repository's main directory, you may run a simulation via

```
julia --project setup_and_run_simulation.jl
```

To make an animation of the simulation, first edit the `visualize_simulation.jl` or `visualize_movie.jl` script to point to the correct `.jld2` output files from the simulation you want to visualize, and then run, e.g.,

```
julia --project visualize_simulation.jl
```

## Cite

You can cite the [GeophysicalFlows.jl](http://github.com/FourierFlows/GeophysicalFlows.jl) package via [zenodo](https://zenodo.org). Please cite as:

> Navid C. Constantinou, Gregory L. Wagner, and co-contributors. (2021). FourierFlows/GeophysicalFlows.jl: GeophysicalFlows v0.11.2  (Version v0.11.2). Zenodo.  [http://doi.org/10.5281/zenodo.1463809](http://doi.org/10.5281/zenodo.1463809)
