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

α = 0.9                       # aspect ratio parameter α = Ly / Lx

ny, Ly  = 256, 2π             # grid resolution and domain length
nx, Lx  = 256, Ly / α         # grid resolution and domain length

## Then we pick the time-stepper parameters
    dt = 1e-2  # timestep
 nsubs = 50    # number of steps between each plot
 
  ν = 1e-4
 k₀ = 1
 tfinal = 0.005 / (ν * k₀^2)
 nsteps = Int(round(tfinal / dt))


# Initial condition (as function)
initial_vorticity(x, y) = sin(α * x) * sin(y) + 0.05 * (cos(α * 8x) + cos(8y))

include("simulation.jl")