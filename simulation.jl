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
# ζ_perturbation = peakedisotropicspectrum(gr, k₀, E₀)
#ζ₀ = @. sin(Y) + ζ_perturbation

X, Y = gridpoints(gr)
ζ₀ = @. initial_vorticity(X, Y) 

TwoDNavierStokes.set_zeta!(prob, ζ₀)

energy_initial = energy(prob)

U = sqrt(2 * energy_initial)

Re = U * L / ν


# ## Diagnostics

"""
    vorticityL4(prob)
Returns the domain-averaged enstrophy, ∫ ζ⁴ dxdy / (Lx Ly), for the solution in `sol`.
"""
@inline function vorticityL4(prob)
  sol, vars, grid = prob.sol, prob.vars, prob.grid
  @. vars.zetah = sol
  ldiv!(vars.zeta, grid.rfftplan, vars.zetah)
  return sum(@. vs.zeta^4) * gr.dx * gr.dy / (gr.Lx * gr.Ly)
end

"""
    palinstrophy(prob)
Returns the domain-averaged palinstrophy, ∫ |∇ζ|² dxdy / (Lx Ly), for the solution in `sol`.
"""
@inline function palinstrophy(prob)
  sol, vars, grid = prob.sol, prob.vars, prob.grid
  palinstrophyh = vars.uh # use vars.uh as scratch variable

  @. palinstrophyh = grid.Krsq * abs2(sol)
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
filepath = "./data/"
filename = joinpath(filepath, "kolmogorovflow.jld2")
filename_diags = joinpath(filepath, "kolmogorovflow_diags.jld2")

filename = FourierFlows.uniquepath(filename)
@info "Output will be saved at $filename."

filename_diags = FourierFlows.uniquepath(filename_diags)
@info "Diagnostics will be saved at $filename_diags."

# Do some basic file management
if isfile(filename); rm(filename); end
if isfile(filename_diags); rm(filename_diags); end

# And then create Output
get_sol(prob) = sol # extracts the Fourier-transformed solution
out = Output(prob, filename, (:zetah, get_sol))
saveproblem(out)


# ## Time-stepping the `Problem` forward

# We time-step the `Problem` forward in time.

saveoutput(out) # save initial condition

file = jldopen(filename, "a+")
FourierFlows.savefield(file, "params/nsubs", nsubs)
FourierFlows.savefield(file, "params/k₀", k₀)
close(file)

@info "Starting simulation..."

startwalltime = time()

for j=0:Int(round(nsteps/nsubs))-1
    
  if j%(1000/nsubs)==0
    TwoDNavierStokes.updatevars!(prob)
    
    cfl = cl.dt * maximum([maximum(vs.u) / gr.dx, maximum(vs.v) / gr.dy])
    
    estimated_remaining_walltime = (time()-startwalltime)/60 / cl.step * (nsteps-cl.step)
    log = @sprintf("step: %04d, t: %d, tν : %.4f, ΔE: %.4f, ΔZL₂: %.4f, ΔZL₄: %.4f, P: %.4f, cfl: %.4f, walltime: %.2f min, estimated remaining walltime: %.2f min",
      cl.step, cl.t, ν*k₀^2*cl.t, E.data[E.i]/E.data[1], Z2.data[Z2.i]/Z2.data[1], Z4.data[Z4.i]/Z4.data[1], P.data[P.i]/P.data[1], cfl, (time()-startwalltime)/60, estimated_remaining_walltime)
    println(log)
  end  

  stepforward!(prob, diags, nsubs)
  dealias!(sol, gr)

  saveoutput(out)
end

@info @sprintf("Simulation finished after %.2f minutes.", (time()-startwalltime)/60)

savediagnostic(E, "energy", filename_diags)
savediagnostic(Z2, "enstrophyL2", filename_diags)
savediagnostic(Z4, "enstrophyL4", filename_diags)
savediagnostic(P, "palinstrophy", filename_diags)

@info "Run visualize_simulation.jl after prescribing the two filenames used for saving output at the top the visualize_simulation.jl script."

@info "Output for this simulation was saved at $filename and $filename_diags."
