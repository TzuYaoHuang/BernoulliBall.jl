using WaterLily,Pathlines,StaticArrays,GLMakie,Plots

# Jet inflow: uniform vertical velocity U inside radius `r` of the domain centreline (dim 2), zero outside
# and zero tangentially — a blower nozzle set into an otherwise closed floor. Passed as the `uBC(i,x,t)`
# function, it also initialises the whole jet column at t=0.
struct Jet{T} <: Function; U::T; xc::T; r::T; end
(bc::Jet)(i,x,t) = i==1 && abs(x[2]-bc.xc)<bc.r ? bc.U : zero(bc.U)

# A ball of diameter L floating in a vertical jet of diameter L, offset `off`·L from the jet centreline.
# The floor (dim-1 low face) is the jet, the sides are slip walls, and the top (dim-1 high face) is a free
# outlet (`exitBC`) — dim 1 must carry the flow for WaterLily's exit convention to apply.
function ball(;L=32,Re=250,U=1,off=0.4,H=10,W=6,T=Float32,mem=Array)
    xc = T(W*L/2)
    center = SA{T}[2L,xc+off*L]
    body = AutoBody((x,t)->√sum(abs2,x-center)-L/2)
    Simulation((H*L,W*L),Jet(T(U),xc,T(1.25L)),L;U,ν=U*L/Re,body,exitBC=true,T,mem)
end

# using CUDA
sim = ball() #;mem=CuArray)
Ni = size(inside(sim.flow.p))
jet = sim.flow.uBC

# pathline tracers (Pathlines.jl): particles seeded across the whole domain, advected with LilyPad.jl's
# 2nd-order departure-point scheme (`Pathlines.update!`, generic over any WaterLily.Simulation) and
# rasterised into a fading, speed-coloured canvas — a numerical dye/smoke visualisation, updated once per
# flow step alongside the pressure field.
particles = Particles(12_000,sim.flow.p;life=UInt(200))
canvas = PathlineCanvas(Ni[1],Ni[2];bgcolor=:black,fadetau=0.8,colormap=:inferno,colorrange=(0,2))

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
GLMakie.contourf!(ax1,p_obs;colormap=:seismic,levels=range(-0.3f0,0.3f0,length=21),extendlow=:auto,extendhigh=:auto)
GLMakie.contour!(ax1,sdf;levels=[0],color=:black,linewidth=2)
GLMakie.image!(ax2,canvas_obs)
GLMakie.contour!(ax2,sdf;levels=[0],color=:white,linewidth=2)

# mark the inflow span on the floor of each panel
for ax ∈ (ax1,ax2)
    GLMakie.lines!(ax,[jet.xc-jet.r,jet.xc+jet.r],[1,1];color=:red,linewidth=6)
    GLMakie.text!(ax,jet.xc,8;text="inflow",align=(:center,:bottom),color=:red,fontsize=13)
end
GLMakie.hidedecorations!(ax1); GLMakie.hidedecorations!(ax2)

# run & visualise
t₀ = sim_time(sim); duration = 80; tstep = 0.2
Cy,t_F = Float64[],Float64[] # lateral force coefficient

GLMakie.record(fig,"ball_on_jet.mp4";framerate=30) do io
    for tᵢ ∈ range(t₀,t₀+duration;step=tstep)
        while sim_time(sim) < tᵢ
            mom_step!(sim.flow,sim.pois)
            dt = Float32(sim.flow.Δt[end-1]) # Δt of the step just taken (Δt[end] is next step's prediction)
            Pathlines.update!(particles,sim)
            fade!(canvas,dt); draw!(canvas,particles.position,particles.position⁰,dt)
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
