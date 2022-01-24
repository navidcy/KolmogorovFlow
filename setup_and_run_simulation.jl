# # 2D decaying turbulence
#
# A simulation of decaying two-dimensional turbulence.

using Pkg; Pkg.instantiate()

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

α = 0.3                       # aspect ratio parameter α = Ly / Lx
γ = 1.0                       # forcing amplitude

add_stochastic_forcing = true
ε = 0.1                       # energy input rate by the forcing

critical_ν(α, γ) = sqrt((γ * sqrt(1 - α^2)) / (sqrt(2) * (1 + α^2)))

ny, Ly  = 64, 2π             # grid resolution and domain length
nx, Lx  = 128, Ly / α

## Then we pick the time-stepper parameters
    dt = 5e-3  # timestep
 nsubs = 50    # number of steps between each plot
 
  ν = critical_ν(α, γ) / 10

 k₀ = 1
 tfinal = 10.0 / (ν * k₀^2)
 nsteps = Int(round(tfinal / dt))

 grid = TwoDGrid(dev, nx, Lx, ny, Ly)
 x, y = gridpoints(grid)

 forcing = @. γ * cos(y)
 forcingh = rfft(forcing)

if add_stochastic_forcing == true
  K = @. sqrt(grid.Krsq)             # a 2D array with the total wavenumber

  forcing_spectrum = 0*grid.Krsq

  forcing_spectrum[grid.Krsq .== 1.0] .= 1
  forcing_spectrum[grid.Krsq .== 4.0] .= 1
  
  forcing_spectrum[grid.Krsq .== 0] .= 0 # ensure forcing has zero domain-average

  ε0 = FourierFlows.parsevalsum(forcing_spectrum .* grid.invKrsq / 2, grid) / (grid.Lx * grid.Ly)
  @. forcing_spectrum *= ε/ε0        # normalize forcing to inject energy at rate ε

  random_uniform = dev==CPU() ? rand : CUDA.rand

  function calcF!(Fh, sol, t, clock, vars, params, grid)
    Fh .= forcingh .+ sqrt.(forcing_spectrum) .* exp.(2π * im * random_uniform(eltype(grid), size(sol))) ./ sqrt(clock.dt)
    return nothing
  end
else
  # Function that constructs the deterministic vorticity forcing in Fourier space
  function calcF!(Fh, sol, t, clock, vars, params, grid)
    @. Fh = forcingh
    return nothing
  end
end

# Initial condition (as function)
initial_vorticity(x, y) = sin(α * x) * sin(y) + 0.05 * (cos(α * 8x) + cos(8y))
initial_vorticity(x, y) = γ / ν * cos(y)

# Setting initial conditions
seed!(1234)

ζ_perturbation = 1 * randn(nx, ny)

ζ₀ = @. initial_vorticity(x, y) + ζ_perturbation

include("simulation.jl")