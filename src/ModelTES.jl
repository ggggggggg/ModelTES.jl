module ModelTES
export
    transitionwidth,
    BiasedTES,
    TESParams,
    getlinearparams,
    NoiseModel,
    TESRecord,
    times,
    rk8,
    IrwinHiltonTES,
    stochastic,
    pulse

using Roots, ForwardDiff, DifferentialEquations
include("rk8.jl")
include("tes_models.jl")

const J_per_eV = 1.602177e-19 #unitless
const kb = 1.38064852e-23 #k boltzmann (J/K)

abstract type AbstractRIT end

# following Irwin-Hilton figure 3
type TESParams{RITType<:AbstractRIT}
    n       ::Float64   # thermal conductance exponent (unitless)
    Tc      ::Float64  # TES critical temperature (K)
    Tbath   ::Float64   # bath temperature (K)

    k       ::Float64   # Pbath pre-factor (W/K^n)
    C       ::Float64   # heat capacity of TES (J/K)

    L       ::Float64  # inductance of SQUID (H)
    Rl      ::Float64   # Thevenin-equivalent resistance Rload = (Rshunt+Rparasitic)(ohms)
                        # note Ibias = V/Rshunt
    Rp      ::Float64   # Rparastic; Rshunt = Rl-Rparasitic (ohms)
    Rn      ::Float64   # normal resistance of TES (ohms)

    RIT     ::RITType   # RIT surface
end


type ShankRIT <: AbstractRIT
    Tw::Float64 # transition width (K)
    A::Float64  # current dependence for R(T,I) (A/K^(3/2))
end
transitionwidth(RIT::ShankRIT)=RIT.Tw
transitionwidth(p::TESParams)=transitionwidth(p.RIT)

"Constructor that fixes `Tw` and `A` for a ShankRIT to have the given `alpha` and `beta`
parameters when biased at resistance `R0`."
function ShankRIT(alpha, beta, n, Tc, Tbath, k, R0, Rn)
    T0 = Tc /(1 + 3*beta/(2*alpha) - 2*(Rn-R0)/(Rn*alpha)*atanh(2*R0/Rn-1))
    I0 = sqrt(k*(T0^n-Tbath^n)/R0)
    Tw = T0*(Rn-R0)/(Rn*log(2)*alpha)
    A = I0*(2*alpha/(3*T0*beta))^(3/2)
    ShankRIT(Tw,A)
end


type BiasedTES{T}
    p::TESParams{T}
    I0::Float64 # intial current for diff eq, aka current through TES (A)
    T0::Float64 # initial temperature for diff equations, aka temperature of TES (K)
    V ::Float64 # Thévinen equivalent voltage V = I0*(p.Rl+p.R0), R0=quiescent resistance
                # also equal to Ibias*Rshunt. Careful! This V is a constant and is NOT
                # the voltage drop across the TES.
end





"For a given T0, the difference R-targetR. Use to solve numerically for T0."
function getR0error(T0, targetR, p::TESParams)
   I0 = sqrt(p.k*(T0^p.n-p.Tbath^p.n)/targetR)
   R(I0,T0,p)-targetR
end


"Find the initial conditions (I0, T0, V) that cause `p` to have resistance `targetR`."
function initialconditions(p::TESParams, targetR)
   T00 = fzero((t)->getR0error(t, targetR,p),p.Tc-10*transitionwidth(p), p.Tc+10*transitionwidth(p))
   I00 = sqrt(p.k*(T00^p.n-p.Tbath^p.n)/targetR)
   R00 = R(I00,T00,p)
   V00 = I00*(p.Rl+R00)
   # now evolve these conditions through integration to really lock them in.
   # shouldn't hard code step size here
   out = pulse(10,1e-1, BiasedTES(p, I00, T00, V00), 0)
   T0 = out.T[end]
   I0 = out.I[end]
   R0 = R(I0,T0,p)
   V = I0*(p.Rl+R0)
   I0,T0,V
end


"Created a biased tes with quiescent state resistance R0"
function BiasedTES{T}(p::TESParams{T}, R0::Float64)
   @assert 0 < R0 < p.Rn
   I0, T0, V = initialconditions(p,R0)
   BiasedTES(p,I0,T0,V)
end

"iv_point(p::TESParams, V, I0, T0)
takes thevinin voltage `V`, and initial current `I0`, and intial temperature `T0`
evolves a pulse for 1 second, and takes the final values
returns I,T,V,R"
function iv_point(p::TESParams, V, I0, T0)
    # solve with an adapative algorithm that is fast for large time steps,
    # ask for very long time steps
    # we probably shouldn't hardcode the time, but 1 second is long for all TESs I know of
    out = pulse(2,1.0, BiasedTES(p, I0, T0, V), 0, method=DifferentialEquations.Rosenbrock23())
    T = out.T[end]
    I = out.I[end]
    R = out.R[end]
    V = I*(p.Rl+R)
    I,T,R,V
end

"iv_curve(p::TESParams, Vs)
takes a sorted array of V thevinin values `Vs`, calculates ivs points by
evolving a pulse for 1 second, and taking the last value
returns Is, Ts, Vs_out, Rs"
function iv_curve(p::TESParams, Vs)
    Is=Vector{Float64}(length(Vs))
    Ts=Vector{Float64}(length(Vs))
    Rs=Vector{Float64}(length(Vs))
    Vs_out=Vector{Float64}(length(Vs))
    @assert issorted(Vs)
    for i in length(Vs):-1:1
            if Vs[i]==0
                I,T,R,V=0.0, p.Tbath, 0.0, 0.0
            elseif i==length(Vs)
                # provide guesses that guaranteed to be resistive
                I,T,R,V = iv_point(p, Vs[i], Vs[i]/p.Rn, p.Tc)
            elseif i<length(Vs)
                #provide last solution as guesses
                # I got a speedup by providing nearby starting points, it was less than a factor of 2
                I,T,R,V = iv_point(p, Vs[i], Is[i+1], Ts[i+1])
            end
            Is[i]=I
            Ts[i]=T
            Rs[i]=R
            Vs_out[i]=V
    end
    Is,Ts,Rs,Vs_out
end

function dT_and_dI_iv_point(p,I,T,R,V)
    bt = BiasedTES(p, I, T, V)
    du = zeros(2)
    u = [T,I]
    bt(0.0, u, du)
    du
end



"Calculate `R0` the quiescent resistance of `tes`."
getR0(tes::BiasedTES) = R(tes.I0, tes.T0, tes.p)


"Calculate thermal conductivity from TES heat capacity to Tbath in small signal limit."
function getG0(tes::BiasedTES)
   p=tes.p
   p.k*p.n*tes.T0^(p.n-1)
end


"Calculate paramaters in Irwin-Hilton table 1."
function getlinearparams(bt::BiasedTES)
   f(x) = R(x[1], x[2], bt.p)
   p = bt.p
   R0 = getR0(bt)
   G0 = getG0(bt)
   drdi,drdt = ForwardDiff.gradient(f, [bt.I0, bt.T0])
   alpha = drdt*bt.T0/R0
   beta = drdi*bt.I0/R0
   PJ = bt.I0^2*R0
   loopgain = PJ*alpha/G0/bt.T0
   tauthermal = p.C/G0
   taucc = tauthermal/(1-loopgain) # constant current time constant
   r = p.Rl/R0
   taueff = (1+beta+r)/(1+beta+r+(1-r)*loopgain) # zero inductance effective thermal time constant
   tauelectrical = p.L/(p.Rl+R0*(1+beta))
   invtau = 1/(2*tauelectrical)+1/(2*taucc)
   a = (1/tauelectrical-1/taucc)^2
   b = -4*(R0/p.L)*loopgain*(2+beta)/tauthermal
   lcritical = -p.L*b/a
   invtaupm = 0.5*sqrt(a+b+0*im) # make it complex, so I can get a complex answer
   tauplus = 1/(invtau+invtaupm)
   tauminus = 1/(invtau-invtaupm)
   c = loopgain*(3+beta-r)+(1+beta+r)
   d = 2*sqrt(loopgain*(2+beta)*(loopgain*(1-r)+(1+beta+r)))
   f = R0*tauthermal/(loopgain-1)^2
   lcritplus = (c+d)*f
   lcritminus = (c-d)*f
   bt.I0, bt.T0, bt.V, p.Rl, p.Tbath, p.Tbath, p.L, R0, G0, p.C, alpha, beta, loopgain,
        tauthermal, taucc, taueff, tauelectrical, tauplus, tauminus, lcritplus, lcritminus, lcritical
end

"Paramters from Irwin-Hilton table one for modeling a linear TES. Defined in Table 1 of Irwin-Hilton chapter."
type IrwinHiltonTES
   I0::Float64
   T0::Float64
   V::Float64
   Rl::Float64
   Tl::Float64 # temperature of the load resistor, usually modeled as =Tbath, but really should come from modeling ep coupling in the load resistor
   Tbath::Float64
   L::Float64
   R0::Float64
   G0::Float64
   C0::Float64
   alpha::Float64
   beta::Float64 #βI
   loopgain::Float64 #ℒI
   tauthermal::Float64
   taucc::Float64 #τI
   taueff::Float64
   tauelectrical::Float64 #τel
   tauplus::Complex{Float64}
   tauminus::Complex{Float64}
   lcritplus::Float64
   lcritminus::Float64
   Lcritical::Float64
end

"Constructor for IrwinHiltonTES model"
function IrwinHiltonTES(tes::BiasedTES)
    linearparams = getlinearparams(tes)
    IrwinHiltonTES(linearparams...)
end

"`Z(tes::IrwinHiltonTES, f)`

Returns the impedance of the `tes` at frequency `f`.
Implements equation 42 of Irwin-Hilton chapter."
function Z(tes::IrwinHiltonTES, f)
  ω=2π*f
  tes.R0*(1+tes.beta) +tes.R0*tes.loopgain*(2+tes.beta)./((1-tes.loopgain)*(1+im*ω*tes.taucc))
end

"`Zcircuit(tes::IrwinHiltonTES, f)`

Returns impedance of complete circuit of `tes` at frequency `f`."
function Zcircuit(tes::IrwinHiltonTES, f)
  ω=2π*f
  tes.R0 + im*ω*tes.L + Z(tes,ω)
end

include("tes_noise.jl")
include("stochastic_integration.jl")

isoverdamped(tes::IrwinHiltonTES) = (isreal(tes.tauplus) && isreal(tes.tauminus) &&
                                    tes.tauplus<tes.tauminus)
isunderdamped(tes::IrwinHiltonTES) = (!isoverdamped(tes) && tes.tauplus!=tes.tauminus)


type TESRecord
    T::Vector{Float64}  # temperature (K)
    I::Vector{Float64}  # current (A)
    R::Vector{Float64}  # TES resistance (Ohm)
    dt::Float64   # seconds between samples (seconds)
end
times(r::TESRecord) = range(0,r.dt,length(r.I))
Base.length(r::TESRecord) = length(r.I)

# I want to re-write this code so that
# A and Tw are the TES params, alpha and beta are derived
# and R0 and be changed to change the bias conditions
# then I0 and T0 and derived from those 3 parametrs

"find bias parameters T0, I0, V for to make the TES have resistane R0 in steady state"
function findbiasparams(tes::TESParams{ShankRIT},R0)

    T0 = Tc /(1 + 3*beta/(2*alpha) - 2*(Rn-R0)/(Rn*alpha)*atanh(2*R0/Rn-1))

    # V = I0*(p.Rl+p.R0)
end


"TES resistance model (Shank et al. 2014)"
R(I, T, RIT::ShankRIT, Tc, Rn) = Rn/2*(1+tanh.((T-Tc+(max.(I,0.0)/RIT.A).^(2/3))/(2*log(2)*RIT.Tw)))
R(I,T, p::TESParams) = R(I,T, p.RIT, p.Tc, p.Rn)


"thermal TES equation
C*dT/dt = I^2*R - k*(T^n-Tbath^n)"
function dT(I, T, k, n, Tbath, C, R)
    # if you run into domain errors, try uncommenting the following line
    # but it probably represents a mistake in your timesteps that should be fixed.
    #T=max(T,0.0) # avoid domain errors when raising T to a power
    Q=-k*(T^n-Tbath^n)+I^2*R
    Q/C

end


"electrical TES equation
L*di/dt = (IBias-I)*Rs+I*Rs-I*Rtes"
function dI(I, T, V, Rl, L, R)
    (V-I*(Rl+R))/L
end

"Calling a BiasedTES gives the dI and dT terms for integration in an in place manner."
function (bt::BiasedTES)(du,u,p_,t) # use DifferentialEquation 4.0+ API
    T,I = u[1],u[2]
    p = bt.p
    r = R(I,T,p)
    # dT(I, T, p.k, p.n, p.Tbath, p.C, r)
    # dI(I,T, bt.V, p.Rl, p.L, r)
    du[1] = dT(I, T, p.k, p.n, p.Tbath, p.C, r)
    du[2] = dI(I,T, bt.V, p.Rl, p.L, r)
end
function (bt::BiasedTES)(t, u, du) # legacy API for rk8
    bt(du,u,nothing,t)
end

function rk8(nsample::Int, dt::Float64, bt::BiasedTES, E::Vector, npresamples::Int=0)
    out = Vector{TESRecord}(length(E))
    for i in 1:length(E)
        out[i] = rk8(nsample, dt, bt, E[i], npresamples)
    end
    out
end

function rk8(nsample::Int, dt::Float64, bt::BiasedTES, E::Number, npresamples::Int=0)
    # Pair of differential equations y' = f(t,y), where y=[T,I]
    p = bt.p
    # Integrate pair of ODEs for all energies EE
    T = Array{Float64}(nsample)
    I = Array{Float64}(nsample)
    T[1:npresamples], I[1:npresamples] = bt.T0, bt.I0 # set T0, I0 for presamples
    y = [bt.T0+E*J_per_eV/p.C, bt.I0]; ys = similar(y); work = Array{Float64}(14)
    T[npresamples+1]=y[1]
    I[npresamples+1]=y[2]
    # npresamples+1 is the point at which initial conditions hold (T differs from T0)
    # npresamples+2 is the first point at which I differs from I0
    for i = npresamples+2:nsample
        rk8!(bt, 0.0, dt, y, ys, work)
        y[:] = ys
        T[i] = y[1]
        I[i] = y[2]
    end
    Rout = [R(I[i],T[i],bt.p) for i=1:length(T)]

    TESRecord(T, I, Rout, dt)
end

# example of using the DifferentialEquations API to solve the relevant equations
function adaptive_solve(bt::BiasedTES, dt::Float64, tspan::Tuple{Float64,Float64}, E::Number, method, abstol, reltol, saveat)
    u0 = [bt.T0+E*ModelTES.J_per_eV/bt.p.C, bt.I0]
    prob = ODEProblem(bt, u0, tspan)
    sol = solve(prob,method,dt=dt,abstol=abstol,reltol=reltol, saveat=saveat, save_everystep=false, dense=false)
end
"pulse(nsample::Int, dt::Float64, bt::BiasedTES, E::Number, npresamples::Int=0; dtsolver=1e-9, method=DifferentialEquations.Tsit5(), abstol=1e-9, reltol=1e-9)"
function pulse(nsample::Int, dt::Float64, bt::BiasedTES, E::Number, npresamples::Int=0; dtsolver=1e-9, method=DifferentialEquations.Tsit5(), abstol=1e-9, reltol=1e-9)
    u0 = [bt.T0+E*ModelTES.J_per_eV/bt.p.C, bt.I0]
    saveat = range(0,dt, nsample-npresamples)
    prob = ODEProblem(bt, u0, (0.0, last(saveat)))
    sol = solve(prob,method,dt=dtsolver,abstol=abstol,reltol=reltol, saveat=saveat, save_everystep=false, dense=false)
    # npresamples+1 is the point at which initial conditions hold (T differs from T0) (sol[1])
    # npresamples+2 is the first point at which I differs from I0
    T = Vector{Float64}(nsample)
    I = Vector{Float64}(nsample)
    T[npresamples+1:end] = sol[1,:]
    I[npresamples+1:end] = sol[2,:]
    T[1:npresamples]=bt.T0
    I[1:npresamples]=bt.I0
    Rout = [R(I[i],T[i],bt.p) for i=1:length(T)]
    TESRecord(T,I, Rout,dt)
end

function pulses(nsample::Int, dt::Float64, bt::BiasedTES, Es::Vector, arrivaltimes::Vector; dtsolver=1e-9, method=DifferentialEquations.Tsit5(), abstol=1e-9, reltol=1e-9)
  u0 = [bt.T0, bt.I0]
  saveat = range(0,dt, nsample)
  prob = ODEProblem(bt, u0, (0.0, last(saveat)))
  Esdict = Dict([(at,E) for (at,E) in zip(arrivaltimes,Es)])
  # this defines a callback that is evaluated when t equals a value in arrival times
  # when evaluated it discontinuously changes the temperature (u[1])
  # the last (true,true) argument has to do with which points are saved
  # it doesn't appear to add extra points beyond saveat, so maybe save_everystep overrides it?
  function cbfun(integrator)
    integrator.u[1]+=Esdict[integrator.t]*ModelTES.J_per_eV/bt.p.C
    # modify the integrator timestep back to dtsolver, to take small steps on the rising edge of the pulse
    integrator.dtpropose=dtsolver # in future use modify_proposed_dt!, see http://docs.juliadiffeq.org/latest/basics/integrator.html#Stepping-Controls-1
  end
  cb = DiscreteCallback((u,t,integrator)->(t in arrivaltimes), cbfun; save_positions=(false,false))
  # tstops is used to make sure the integrator checks each time in arrivaltimes
  sol = solve(prob,method,dt=dtsolver,abstol=abstol,reltol=reltol, saveat=saveat, save_everystep=false, dense=false,callback=cb, tstops=arrivaltimes)

  T = sol[1,:]
  I = sol[2,:]
  Rout = [R(I[i],T[i],bt.p) for i=1:length(T)]
  TESRecord(T,I, Rout,dt)
end

end # module
