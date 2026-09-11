using WaterLily,BiotSavartBCs,StaticArrays,GLMakie,Plots,Adapt

# Jet inflow: uniform axial velocity U inside radius `r` of the domain's (dim2,dim3) centreline,
# zero outside and zero tangentially — a circular blower nozzle set into an otherwise closed floor.
# Passed as the `uBC(i,x,t)` function, it also initialises the whole jet column at t=0.
struct Jet{T} <: Function; U::T; xc::T; r::T; end
(bc::Jet)(i,x,t) = i==1 && √((x[2]-bc.xc)^2+(x[3]-bc.xc)^2)<bc.r ? bc.U : zero(bc.U)
# Jet's 3 fields share a single type parameter T, which breaks Adapt's generic closure adaptor
# (it assumes one type parameter per captured field). Since the fields are plain scalars that
# never need device conversion, adapting is just the identity.
Adapt.adapt_structure(to,j::Jet) = j

# A ball of diameter L floating in a vertical jet of diameter L, offset `off`·L from the jet axis
# along dim2 only (the domain's (dim2,dim3) cross-section is square, so the jet nozzle sits at its
# centre). The floor (dim-1 low face) is the jet inflow; the top and sides are open far-field
# boundaries handled by BiotSavartBCs.jl (velocity there is set from the interior vorticity via the
# Biot-Savart integral, rather than a slip wall or convective exit), so only the floor (`-1`) is
# excluded via `nonbiotfaces`.
function ball(;D=2^5,Re=10^5,U=1,off=0.3,H=4,W=6,h=2,r=1.5,T=Float32,mem=Array)
    R = T(D/2)
    xc = T(W*R)
    center = SA{T}[h*D,xc+off*R,xc]
    body = AutoBody((x,t)->√sum(abs2,x-center)-R)
    BiotSimulation((H*D,W*D,W*D),Jet(T(U),xc,T(r*R)),D;U,ν=U*D/Re,body,T,mem,nonbiotfaces=(-1,))
end

using CUDA
mem = CUDA.functional() ? CuArray : Array
sim = ball(;mem)
t₀ = sim_time(sim); duration = 100; tstep = 0.2

# λ2 criterion (Jeong & Hussain 1995): iso-surfaces of log10(-λ2) pick out vortex cores, unlike
# vorticity magnitude which also lights up shear layers. Same recipe as WaterLily-Examples'
# ThreeD_TaylorGreenVortex.jl — computed into sim.flow.σ, then copied to the CPU buffer viz! renders.
# `viz!` calls this `f` once per rendered frame (after stepping to that frame's time), with a fully
# projected, divergence-free flow field — unlike `udf`, which fires mid-step on unprojected data — so
# it also doubles as the hook for sampling the lateral force coefficient once per frame.
Cy,t_F = Float64[],Float64[] # lateral force coefficient, dim2 being the ball's offset axis
function λ₂!(arr,sim)
    a = sim.flow.σ
    @inside a[I] = log10(max(1e-6,-WaterLily.λ₂(I,sim.flow.u)*sim.L/sim.U))
    copyto!(arr,a[inside(a)])
    push!(Cy,-2WaterLily.pressure_force(sim)[2]/sim.L); push!(t_F,sim_time(sim))
end
# no `sym` here (unlike TGV/jelly): the ball's off-axis offset breaks the domain's mirror symmetry
viz!(sim;f=λ₂!,duration,step=tstep,algorithm=:absorption,colormap=:Reds,video="ball_on_jet_3d.mp4")

# lateral force coefficient vs time, sampled once per rendered frame (coarser than the 2D script,
# which samples every internal solver step) — see TwoD_BallOnJet.jl for how to read this curve.
Plots.plot(t_F,Cy,xlabel="tU/L",ylabel="Cy",legend=false)
savefig("restoring_force_3d.png")
