using WaterLily,BiotSavartBCs,Pathlines,StaticArrays,GLMakie,Plots,Adapt

# Jet inflow: uniform vertical velocity U inside radius `r` of the domain centreline (dim 2), zero outside
# and zero tangentially — a blower nozzle set into an otherwise closed floor. Passed as the `uBC(i,x,t)`
# function, it also initialises the whole jet column at t=0.
struct Jet{T} <: Function; U::T; xc::T; r::T; end
(bc::Jet)(i,x,t) = i==1 && abs(x[2]-bc.xc)<bc.r ? bc.U : zero(bc.U)
# Jet's 3 fields share a single type parameter T, which breaks Adapt's generic closure adaptor
# (it assumes one type parameter per captured field). Since the fields are plain scalars that
# never need device conversion, adapting is just the identity.
Adapt.adapt_structure(to,j::Jet) = j

# A ball of diameter L floating in a vertical jet of diameter L, offset `off`·L from the jet centreline.
# The floor (dim-1 low face) is the jet inflow; the top and sides are open far-field boundaries handled
# by BiotSavartBCs.jl (velocity there is set from the interior vorticity via the Biot-Savart integral,
# rather than a slip wall or convective exit), so only the floor (`-1`) is excluded via `nonbiotfaces`.
function ball(;D=2^5,Re=10^5,U=1,off=0.3,H=4,W=6,h=2,r=1.5,T=Float32,mem=Array)
    R = T(D/2)
    xc = T(W*R)
    center = SA{T}[h*D,xc+off*R]
    # body = WaterLily.NoBody() 
    body = AutoBody((x,t)->√sum(abs2,x-center)-R)
    BiotSimulation((H*D,W*D),Jet(T(U),xc,T(r*R)),D;U,ν=U*D/Re,body,T,mem,nonbiotfaces=(-1,))
end

using CUDA
mem = CUDA.functional() ? CuArray : Array
sim = ball(;mem)
Ni = size(inside(sim.flow.p))
jet = sim.flow.uBC

# pathline tracers (Pathlines.jl): particles seeded across the whole domain, advected with LilyPad.jl's
# 2nd-order departure-point scheme (`Pathlines.update!`, generic over any WaterLily.Simulation) and
# rasterised into a fading, speed-coloured canvas — a numerical dye/smoke visualisation, updated once per
# flow step alongside the pressure field.
particles = Particles(16_000,sim.flow.p;life=UInt(200),mem)
canvas = PathlineCanvas(Ni[1],Ni[2];bgcolor=:white,fadetau=1.2,colormap=:inferno,colorrange=(0,.9))

# ball outline: computed once since the body is static. Raw pressure inside the immersed body is not
# physically meaningful, so it's masked out (NaN) wherever the sdf is negative.
WaterLily.measure_sdf!(sim.flow.σ,sim.body,WaterLily.time(sim.flow))
sdf = permutedims(Array(sim.flow.σ[inside(sim.flow.σ)]))
masked_p() = ifelse.(sdf.<0,NaN32,permutedims(Array(sim.flow.p[inside(sim.flow.p)])))

# figure: pressure field + ball on the left, pathline canvas + ball on the right. Both share the sim's
# (dim1,dim2)=(z,x) grid, transposed with `permutedims` so the vertical jet plots vertically (Makie plots
# a matrix's 1st dimension along x, 2nd along y — dim 1 here is the flow direction, not x).
p_obs = Observable(masked_p())
canvas_obs = Observable(permutedims(canvas.canvas))

# GLMakie and Plots both export plotting verbs (contourf!, contour!, image!, ...), so qualify the Makie
# calls explicitly and keep Plots only for the final force plot below
fig = GLMakie.Figure(size=(1000,650))
ax1 = GLMakie.Axis(fig[1,1],aspect=GLMakie.DataAspect(),title="pressure")
ax2 = GLMakie.Axis(fig[1,2],aspect=GLMakie.DataAspect(),title="pathlines")
# alpha<1 softens the colormap so the red/blue extremes aren't fully saturated, toning down the contrast
GLMakie.contourf!(ax1,p_obs;colormap=(:seismic,0.6),levels=range(-0.3f0,0.3f0,length=21),extendlow=:auto,extendhigh=:auto)
GLMakie.contour!(ax1,sdf;levels=[0],color=:black,linewidth=2)
# canvas_obs is rendered at its own fixed pixel resolution, independent of the sim's grid size — so it
# must be stretched explicitly onto the grid's coordinate range (0..Ni[2] × 0..Ni[1], matching sdf below)
# or it plots at pixel-count scale instead, making the body/inlet markers look wildly mis-sized against it
GLMakie.image!(ax2,0..Ni[2],0..Ni[1],canvas_obs)
# fill the body solid so it reads as an object against the pathlines, not just a thin outline
GLMakie.contourf!(ax2,sdf;levels=[-1f4,0],colormap=[:gray50])
GLMakie.contour!(ax2,sdf;levels=[0],color=:black,linewidth=2)

# mark the inflow span on the floor of each panel
for ax ∈ (ax1,ax2)
    GLMakie.lines!(ax,[jet.xc-jet.r,jet.xc+jet.r],[1,1];color=:red,linewidth=6)
    GLMakie.text!(ax,jet.xc,8;text="inflow",align=(:center,:bottom),color=:red,fontsize=13)
end
GLMakie.hidedecorations!(ax1); GLMakie.hidedecorations!(ax2)

# run & visualise
t₀ = sim_time(sim); duration = 160; tstep = 0.2
Cy,t_F = Float64[],Float64[] # lateral force coefficient

GLMakie.record(fig,"ball_on_jet.mp4";framerate=30) do io
    for tᵢ ∈ range(t₀,t₀+duration;step=tstep)
        while sim_time(sim) < tᵢ
            mom_step!(sim.flow,sim.pois)
            dt = Float32(sim.flow.Δt[end-1]) # Δt of the step just taken (Δt[end] is next step's prediction)
            Pathlines.update!(particles,sim)
            fade!(canvas,dt); draw!(canvas,Array(particles.position),Array(particles.position⁰),dt)
            push!(Cy,-2WaterLily.pressure_force(sim)[2]/sim.L); push!(t_F,sim_time(sim))
        end
        isempty(Cy) || println("tU/L=",round(tᵢ,digits=3),"  Cy=",round(Cy[end],digits=3))
        p_obs[] = masked_p()
        canvas_obs[] = permutedims(canvas.canvas)
        recordframe!(io)
    end
end

# lateral force coefficient: sharply negative (pulling the off-centre ball back toward the jet axis) as
# the Bernoulli suction hits on startup, then oscillates with the shed wake — its long-time mean is small
# and sign-sensitive to `off` and Re, rather than a settled restoring constant.
Plots.plot(t_F,Cy,xlabel="tU/L",ylabel="Cy",legend=false)
savefig("restoring_force.png")
