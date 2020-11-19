# # 2D decaying turbulence
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


# ## Numerical, domain, and simulation parameters
#
# First, we pick some numerical and physical parameters for our model.

n, L  = 128, 2π             # grid resolution and domain length

## Then we pick the time-stepper parameters
    dt = 5e-2  # timestep
 nsubs = 20    # number of steps between each plot
 
  ν = 1e-4
 k₀ = 1
 tfinal = 0.02 / (ν * k₀^2)
 
 nsteps = Int(round(tfinal / dt))


# ## Problem setup
# We initialize a `Problem` by providing a set of keyword arguments. The
# `stepper` keyword defines the time-stepper to be used.
prob = TwoDNavierStokes.Problem(dev; nx=n, Lx=L, ny=n, Ly=L, ν=ν, dt=dt, stepper="ETDRK4")

# Next we define some shortcuts for convenience.
sol, cl, vs, gr = prob.sol, prob.clock, prob.vars, prob.grid
x, y = gr.x, gr.y


# ## Setting initial conditions

# Our initial condition closely tries to reproduce the initial condition used
# in the paper by McWilliams (_JFM_, 1984)
seed!(1234)
E₀ = 0.00025 # energy of sin(y) is 0.25

# ζ_perturbation = peakedisotropicspectrum(gr, k₀, E₀, mask=prob.timestepper.filter)
ζ_perturbation = peakedisotropicspectrum(gr, k₀, E₀)
X, Y = gridpoints(gr)
# ζ₀ = @. sin(Y) + ζ_perturbation
ζ₀ = @. sin(X) * sin(Y) + 0.005 * cos(X)

TwoDNavierStokes.set_zeta!(prob, ζ₀)

energy_initial = energy(prob)
U = sqrt(2*energy_initial)
Re = U * L / ν


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


# ## Output

# We choose folder for outputing `.jld2` files and snapshots (`.png` files).
filepath = "."
plotpath = "./plots_barotropic_gammaplane"
plotname = "snapshots"
filename = joinpath(filepath, "kolmogorovflow.jld2")
filename_diags = joinpath(filepath, "kolmogorovflow_diags.jld2")

filename = FourierFlows.uniquepath(filename)
@info "Output will be saved at $filename."

filename_diags = FourierFlows.uniquepath(filename_diags)
@info "Diagnostics will be saved at$filename_diags."

# Do some basic file management
if isfile(filename); rm(filename); end
if isfile(filename_diags); rm(filename_diags); end
if !isdir(plotpath); mkdir(plotpath); end

# And then create Output
get_sol(prob) = sol # extracts the Fourier-transformed solution
out = Output(prob, filename, (:zetah, get_sol))
saveproblem(out)

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
              legend = :topright,
           linewidth = 2,
               alpha = 0.7,
              xlabel = "ν k₀² t",
               xlims = (0, 1.01 * ν * k₀^2 * tfinal),
               ylims = (0, 3))

l = @layout Plots.grid(1, 2)
p = plot(p1, p2, layout = l, size = (900, 400))


# ## Time-stepping the `Problem` forward

# We time-step the `Problem` forward in time.

startwalltime = time()

saveoutput(out) # save initial condition

@info "Starting simulation..."

for j=0:Int(nsteps/nsubs)-1
  
  cfl = cl.dt * maximum([maximum(vs.u) / gr.dx, maximum(vs.v) / gr.dy])
  
  if j%(1000/nsubs)==0
    log = @sprintf("step: %04d, t: %d, tν : %.4f, ΔE: %.4f, ΔZL₂: %.4f, ΔZL₄: %.4f, P: %.4f, cfl: %.4f, walltime: %.2f min",
      cl.step, cl.t, ν*k₀^2*cl.t, E.data[E.i]/E.data[1], Z2.data[Z2.i]/Z2.data[1], Z4.data[Z4.i]/Z4.data[1], P.data[P.i]/P.data[1], cfl, (time()-startwalltime)/60)
    println(log)
  end  

  stepforward!(prob, diags, nsubs)

  saveoutput(out)
end

@info @sprintf("Simulation finished after %.2f minutes.", (time()-startwalltime)/60)

savediagnostic(E, "energy", filename_diags)
savediagnostic(Z2, "enstrophyL2", filename_diags)
savediagnostic(Z4, "enstrophyL4", filename_diags)
savediagnostic(P, "palinstrophy", filename_diags)

@info "Run visualize_simulation.jl after prescribing the two filenames above in the top of the script."