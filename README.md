# KolmogorovFlow

Study of the stability of Kolmogorov flow. The code is written in [Julia](https://julialang.org).

## Installation

To install the code first you need to [install Julia](https://julialang.org/downloads/). Then glone the repository, e.g.,

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

To create an animation of the simulation, first edit the `visualize_simulation.jl` or `visualize_movie.jl` scripts to point to the correct `.jld2` output file from the simulation you want to visualize, and then run, e.g.,


```
julia --project visualize_simulation.jl
```
