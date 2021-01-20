# KolmogorovFlow

Study of the stability of Kolmogorov flow. The code is written in [Julia](https://julialang.org).

## Installation

First you need to [install Julia](https://julialang.org/downloads/). Then clone the repository, e.g.,

```
git clone https://github.com/navidcy/KolmogorovFlow.git
```

Enter the directory you've cloned the repository and create two directories where simulation output and animations will be saved at. The code will expect these directories to be called `data` and `movies`, i.e.,

```
mkdir movies data
```

Then, after you edit `setup_and_run_simulation.jl` file with your parameter values, and while still inside the repository's main directory, you may run a simulation via

```
julia --project setup_and_run_simulation.jl
```

To make an animation of the simulation, first edit the `visualize_simulation.jl` or `visualize_movie.jl` script to point to the correct `.jld2` output files from the simulation you want to visualize, and then run, e.g.,

```
julia --project visualize_simulation.jl
```
