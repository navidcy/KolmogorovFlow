# # 2D decaying turbulence
#
# A simulation of decaying two-dimensional turbulence.

using FourierFlows, Printf, Random, JLD2
 
using Random: seed!
using FFTW: rfft, irfft
using LinearAlgebra: ldiv!
using FourierFlows: parsevalsum

import GeophysicalFlows.TwoDNavierStokes
import GeophysicalFlows.TwoDNavierStokes: energy, enstrophy
import GeophysicalFlows: peakedisotropicspectrum


# ## Choosing a device: CPU or GPU

dev = CPU()     # Device (CPU/GPU)

# ## Numerical, domain, and simulation parameters
#
# First, we pick some numerical and physical parameters for our model.

α = 1.2                       # aspect ratio parameter α = Ly / Lx

ny, Ly  = 128, 2π             # grid resolution and domain length
nx, Lx  = 126, Ly / α

## Then we pick the time-stepper parameters
    dt = 1e-3  # timestep
 nsubs = 1    # number of steps between each plot
 
  ν = 2e-1
 k₀ = 1
 tfinal = 0.2 / (ν * k₀^2)
 nsteps = Int(round(tfinal / dt))

 
 grid = TwoDGrid(dev, nx, Lx, ny, Ly)
 x, y = gridpoints(grid)

 γ = 1.0
 forcing = @. γ * cos(y)
 forcingh = rfft(forcing)

 # Function that constructs the vorticity forcing in Fourier space
 function calcF!(Fh, sol, t, clock, vars, params, grid)
   @. Fh = forcingh
   return nothing
 end


# Initial condition (as function)
initial_vorticity(x, y) = sin(α * x) * sin(y) + 0.05 * (cos(α * 8x) + cos(8y))
initial_vorticity(x, y) = γ / ν * cos(y)

# Setting initial conditions
seed!(1234)

ζ_perturbation = 1 * randn(nx, ny)

ζ₀ = @. initial_vorticity(x, y) + ζ_perturbation

include("simulation.jl")