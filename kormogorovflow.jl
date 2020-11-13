# # 2D decaying turbulence
#
#md # This example can be run online via [![](https://mybinder.org/badge_logo.svg)](@__BINDER_ROOT_URL__/generated/twodnavierstokes_decaying.ipynb).
#md # Also, it can be viewed as a Jupyter notebook via [![](https://img.shields.io/badge/show-nbviewer-579ACA.svg)](@__NBVIEWER_ROOT_URL__/generated/twodnavierstokes_decaying.ipynb).
#
# A simulation of decaying two-dimensional turbulence.

using FourierFlows, Printf, Random, Plots
 
using Random: seed!
using FFTW: rfft, irfft
using LinearAlgebra: ldiv!
using FourierFlows: parsevalsum

import GeophysicalFlows.TwoDNavierStokes
import GeophysicalFlows.TwoDNavierStokes: energy, enstrophy
import GeophysicalFlows: peakedisotropicspectrum


# ## Choosing a device: CPU or GPU

dev = CPU()     # Device (CPU/GPU)
nothing # hide


# ## Numerical, domain, and simulation parameters
#
# First, we pick some numerical and physical parameters for our model.

n, L  = 512, 2π             # grid resolution and domain length
nothing # hide

## Then we pick the time-stepper parameters
    dt = 1e-2  # timestep
 nsubs = 80    # number of steps between each plot
 
  ν = 3e-5
 k₀ = 10
 tfinal = 2 / (ν * k₀^2)
 
 nsteps = Int(round(tfinal / dt))
 
nothing # hide


# ## Problem setup
# We initialize a `Problem` by providing a set of keyword arguments. The
# `stepper` keyword defines the time-stepper to be used.
prob = TwoDNavierStokes.Problem(dev; nx=n, Lx=L, ny=n, Ly=L, ν=ν, dt=dt, stepper="ETDRK4")
nothing # hide

# Next we define some shortcuts for convenience.
sol, cl, vs, gr = prob.sol, prob.clock, prob.vars, prob.grid
x, y = gr.x, gr.y
nothing # hide


# ## Setting initial conditions

# Our initial condition closely tries to reproduce the initial condition used
# in the paper by McWilliams (_JFM_, 1984)
seed!(1234)
E₀ = 0.00025 # energy of sin(y) is 0.25

# ζ_perturbation = peakedisotropicspectrum(gr, k₀, E₀, mask=prob.timestepper.filter)
ζ_perturbation = peakedisotropicspectrum(gr, k₀, E₀)
X, Y = gridpoints(gr)
# ζ₀ = @. sin(Y) + ζ_perturbation
ζ₀ = @. sin(X) * sin(Y) + 0.1 * sin(2 * X) * sin(2 * Y)

TwoDNavierStokes.set_zeta!(prob, ζ₀)

energy_initial = energy(prob)
U = sqrt(2*energy_initial)
Re = U * L / ν

nothing # hide

# Let's plot the initial vorticity field:
p = heatmap(x, y, vs.zeta',
         aspectratio = 1,
              c = :balance,
           clim = (-2, 2),
          xlims = (-L/2, L/2),
          ylims = (-L/2, L/2),
         xticks = -3:3,
         yticks = -3:3,
         xlabel = "x",
         ylabel = "y",
          title = "initial vorticity",
     framestyle = :box)


# ## Diagnostics

"""
    vorticityL4(prob)
Returns the domain-averaged enstrophy, ∫ ½ ζ⁴ dxdy / (Lx Ly), for the solution in `sol`.
"""
@inline function vorticityL4(prob)
  sol, vars, grid = prob.sol, prob.vars, prob.grid
  @. vars.zetah = sol
  ldiv!(vars.zeta, grid.rfftplan, vars.zetah)
  return sum(@. 0.5 * vs.zeta^4) * gr.dx * gr.dy / (gr.Lx * gr.Ly)
end

"""
    palinstrophy(prob)
Returns the domain-averaged palinstrophy, ∫ ½ |∇ζ|² dxdy / (Lx Ly), for the solution in `sol`.
"""
@inline function palinstrophy(prob)
  sol, vars, grid = prob.sol, prob.vars, prob.grid
  palinstrophyh = vars.uh # use vars.uh as scratch variable

  @. palinstrophyh = 1 / 2 * grid.Krsq * abs2(sol)
  return 1 / (grid.Lx * grid.Ly) * parsevalsum(palinstrophyh, grid)
end


# Create Diagnostics -- `energy` and `enstrophy` functions are imported at the top.
E = Diagnostic(energy, prob; nsteps=nsteps)
Z2 = Diagnostic(enstrophy, prob; nsteps=nsteps)
Z4 = Diagnostic(vorticityL4, prob; nsteps=nsteps)
P = Diagnostic(palinstrophy, prob; nsteps=nsteps)
diags = [E, Z2, Z4, P] # A list of Diagnostics types passed to "stepforward!" will  be updated every timestep.
nothing # hide


# ## Output

# We choose folder for outputing `.jld2` files and snapshots (`.png` files).
filepath = "."
plotpath = "./plots_decayingTwoDNavierStokes"
plotname = "snapshots"
filename = joinpath(filepath, "decayingTwoDNavierStokes.jld2")
nothing # hide

# Do some basic file management
if isfile(filename); rm(filename); end
if !isdir(plotpath); mkdir(plotpath); end
nothing # hide

# And then create Output
get_sol(prob) = Array(prob.sol) # extracts the Fourier-transformed solution
out = Output(prob, filename, (:sol, get_sol))
saveproblem(out)
nothing # hide


# ## Visualizing the simulation

# We initialize a plot with the vorticity field and the time-series of
# energy and enstrophy diagnostics.

p1 = heatmap(x, y, vs.zeta',
         aspectratio = 1,
                   c = :balance,
                clim = (-2, 2),
               xlims = (-L/2, L/2),
               ylims = (-L/2, L/2),
              xticks = -3:3,
              yticks = -3:3,
              xlabel = "x",
              ylabel = "y",
               title = "vorticity, t="*@sprintf("%.2f", cl.t),
          framestyle = :box)

p2 = plot(4, # this means "a plot with two series"
               label = ["|u|_{L2}/|ζ|_{L2}(t=0)"  "|ζ|_{L2}/|ζ|_{L2}(t=0)" "|ζ|_{L4}/|ζ|_{L4}(t=0)" "|∇ζ|_{L2}/|∇ζ|_{L2}(t=0)"],
              legend = :bottomleft,
           linewidth = 2,
               alpha = 0.7,
              xlabel = "ν k₀² t",
               xlims = (0, 1.01 * ν * k₀^2 * tfinal),
               ylims = (0, 2))

l = @layout grid(1, 2)
p = plot(p1, p2, layout = l, size = (900, 400))


# ## Time-stepping the `Problem` forward

# We time-step the `Problem` forward in time.

startwalltime = time()

anim = @animate for j = 0:Int(round(nsteps/nsubs))
  
  cfl = cl.dt * maximum([maximum(vs.u) / gr.dx, maximum(vs.v) / gr.dy])
  
  log = @sprintf("step: %04d, t: %d, tν : %.4f, ΔE: %.4f, ΔZ2: %.4f, ΔZ4: %.4f, cfl: %.4f, walltime: %.2f min",
      cl.step, cl.t, ν*k₀^2*cl.t, E.data[E.i]/E.data[1], Z2.data[Z2.i]/Z2.data[1], Z4.data[Z4.i]/Z4.data[1], cfl, (time()-startwalltime)/60)
  
  if j%(1000/nsubs)==0; println(log) end  
  
  p[1][1][:z] = vs.zeta'
  p[1][:title] = "vorticity, ν k₀² t = "*@sprintf("%.2f", ν*k₀^2*cl.t)
  p[2][:title] = "Re = "*@sprintf("%.2f", Re)
  push!(p[2][1], ν * k₀^2 * E.t[E.i], (E.data[E.i]/E.data[1])^(1/2))
  push!(p[2][2], ν * k₀^2 * Z2.t[Z2.i], (Z2.data[Z2.i]/Z2.data[1])^(1/2))
  push!(p[2][3], ν * k₀^2 * Z4.t[Z4.i], (Z4.data[Z4.i]/Z4.data[1])^(1/4))
  push!(p[2][4], ν * k₀^2 * P.t[P.i], (P.data[P.i]/P.data[1])^(1/2))

  stepforward!(prob, diags, nsubs)
  TwoDNavierStokes.updatevars!(prob)  
  
end

gif(anim, "kolmogorov.gif", fps=18)
mp4(anim, "kolmogorov.mp4", fps=18)
