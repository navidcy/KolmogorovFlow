using JLD2
using FourierFlows
using Plots
using Printf
using LinearAlgebra: ldiv!

filename = "./kolmogorovflow_10.jld2"
filename_diags = "./kolmogorovflow_diags_6.jld2"

withoutgif(path) = (length(path)>3 && path[end-3:end] == ".gif") ? path[1:end-4] : path

"""
    uniquepath(path)

Returns `path` with a number appended if `isfile(path)`, incremented until `path` does not exist.
"""
function uniquegifpath(path)
  n = 1
  if isfile(path)
    path = withoutgif(path) * "_$n.gif"
  end
  while isfile(path)
    n += 1
    path = withoutgif(path)[1:end-length("_$(n-1)")] * "_$n.gif"
  end
  path
end

moviegif_filename = "./kolmogorov.gif"
moviegif_filename = uniquegifpath(moviegif_filename)
moviemp4_filename = moviegif_filename[1:end-4]*".mp4"

# ## Visualizing the simulation

function plot_output(x, y, ζ, ψ, t, k₀, ν, t_final)
    tν = ν * k₀^2 * t
    
    p_ζ = heatmap(x, y, ζ',
             aspectratio = 1,
                       c = :balance,
                    clim = (-2, 2),
                   xlims = (x[1], x[end]),
                   ylims = (y[1], y[end]),
                  xticks = -3:3,
                  yticks = -3:3,
                  xlabel = "x",
                  ylabel = "y",
                   title = "vorticity, tν = "*@sprintf("%.2f", tν),
              framestyle = :box)

    p_ψ = contourf(x, y, ψ',
               aspectratio = 1,
                    levels = 11,
                     xlims = (x[1], x[end]),
                     ylims = (y[1], y[end]),
                    xticks = -3:3,
                    yticks = -3:3,
                    xlabel = "x",
                    ylabel = "y",
                     title = "streamfunction",
                framestyle = :box)

    p_diags = plot(5, # this means "a plot with two series"
                   label = ["|u|₂ / |u|₂(t=0)"  "|ζ|₂ / |ζ|₂(t=0)" "|ζ|₄ / |ζ|₄(t=0)" "|∇ζ|₂ / |∇ζ|₂(t=0)"],
                  legend = :topright,
               linewidth = 2,
                   alpha = 0.7, 
                  xlabel = "ν k₀² t",
                   xlims = (0, 1.01 * ν * k₀^2 * t_final),
                   # ylims = (0, 3),
                  yscale = :log10,
                   )

    l = @layout [ Plots.grid(1, 2)
                  c{0.4h} ]
             
    p = plot(p_ζ, p_ψ, p_diags, layout = l, size = (900, 700))

    return p
end

diags = jldopen(filename_diags)


# now let's make a movie of the flow fields

file = jldopen(filename)

nx, ny = file["grid/nx"], file["grid/ny"]
Lx, Ly = file["grid/Lx"], file["grid/Ly"]

global ν = file["params/ν"]
global k₀ = file["params/k₀"]
global nsubs = file["params/nsubs"]

grid = TwoDGrid(nx, Lx, ny, Ly)

x, y = grid.x, grid.y

iterations = parse.(Int, keys(file["snapshots/t"]))
final_iteration = iterations[end]

t_final = file["snapshots/t/$final_iteration"]

ζ = zeros(nx, ny)
ψ = zeros(nx, ny)
ψh = zeros(Complex{Float64}, (grid.nkr, grid.nl))

ζh_initial = file["snapshots/zetah/0"]

energyh = @. 1 / 2 * grid.invKrsq * abs2(ζh_initial)
energy_initial = 1 / (grid.Lx * grid.Ly) * FourierFlows.parsevalsum(energyh, grid)
U = sqrt(2*energy_initial)

global Re = U * grid.Lx / ν

global p = plot_output(x, y, ζ, ψ, 0.0, k₀, ν, t_final)

anim = @animate for (i, iteration) in enumerate(iterations)
  t = file["snapshots/t/$iteration"]
  tν = ν * k₀^2 * t
  
  ζh = file["snapshots/zetah/$iteration"]
  
  @. ψh = - grid.invKrsq * ζh
  ldiv!(ζ, grid.rfftplan, ζh)
  ldiv!(ψ, grid.rfftplan, ψh)
  
  p[1][1][:z] = ζ'
  p[1][:title] = "vorticity, tν = "*@sprintf("%.4f", tν)
  p[2][1][:z] = ψ'
  p[2][:title] = "streamfunction"
  p[3][:title] = "Re = "*@sprintf("%.2f", Re)
  
  t = diags["diags/energy/t"][(i-1)*nsubs+1]
  tν = ν * k₀^2 * t
  
  ΔE = diags["diags/energy/data"][(i-1)*nsubs+1] / diags["diags/energy/data"][1]
  ΔΖ₂ = diags["diags/enstrophyL2/data"][(i-1)*nsubs+1] / diags["diags/enstrophyL2/data"][1]
  ΔΖ₄ = diags["diags/enstrophyL4/data"][(i-1)*nsubs+1] / diags["diags/enstrophyL4/data"][1]
  ΔP = diags["diags/palinstrophy/data"][(i-1)*nsubs+1] / diags["diags/palinstrophy/data"][1]
  psixy00 = diags["diags/psixy00/data"][(i-1)*nsubs+1]
  
  push!(p[3][1], tν, (ΔE)^(1/2))
  push!(p[3][2], tν, (ΔΖ₂)^(1/2))
  push!(p[3][3], tν, (ΔΖ₄)^(1/4))
  push!(p[3][4], tν, (ΔP)^(1/2))
  push!(p[3][5], tν, psixy00)
end

gif(anim, moviegif_filename, fps=14)
mp4(anim, moviemp4_filename, fps=14)
