module ModelTES
export
    BiasedTES,
    TESParams,
    ShankRIT,
    pulses

using Roots, ForwardDiff, DifferentialEquations, Unitful, Base.Test

abstract type AbstractRIT end

# following Irwin-Hilton figure 3
type TESParams{RITType<:AbstractRIT}
    n       ::Float64   # thermal conductance exponent (unitless)
    Tbath   ::typeof(1.0u"mK")  # bath temperature

    G       ::typeof(1.0u"pW/K")  # thermal conductivity G = n*k*(T^n-1)
    C       ::typeof(1.0u"pJ/K")   # heat capacity of TES

    L       ::typeof(1.0u"nH") # inductance of SQUID (H)
    Rl      ::typeof(1.0u"mΩ")  # Thevenin-equivalent resistance Rload = (Rshunt+Rparasitic)(ohms)
                        # note Ibias = V/Rshunt
    Rp      ::typeof(1.0u"mΩ")    # Rparastic; Rshunt = Rl-Rparasitic (ohms)

    RIT     ::RITType   # RIT surface
end
transitiontemperature(p::TESParams) = transitiontemperature(p.RIT)

"thermal_transport_prefactor(p::TESParams)
Return κ of P=κ(T^n-Tbath^n) calculated from p.G and transitiontemperature(G)"
function thermal_transport_prefactor(p::TESParams)
    nn = round(Int,n*100)//100 # make a rational number with a reasonable number of digits
    # since Unitful only supports rational exponents
    uconvert(u"nW/K",p.G)/n/uconvert(u"K",transitiontemperature(p))^(nn-1)
end


type ConstantRIT <: AbstractRIT
    Rn::typeof(1.0u"mΩ")
end
R(I,T,RIT::ConstantRIT) = RIT.Rn
type ShankRIT <: AbstractRIT
    Tc::typeof(1.0u"mK")
    Rn::typeof(1.0u"mΩ")
    Tw::typeof(1.0u"mK")  # transition width
    A::typeof(1.0u"A/K^(3/2)")  # current dependence for R(T,I) (A/K^(3/2))
end
"R(I, T, RIT::ShankRIT)
Return resistance as a function of current `I` and temperature `T` using the model `RIT`.
TES resistance model (Shank et al. 2014)"
R(I, T, RIT::ShankRIT) = RIT.Rn/2*(1+tanh.(Float64.((T-RIT.Tc+(max.(I,0.0u"mA")/RIT.A).^(2/3))/(2*log(2)*RIT.Tw))))
transitionwidth(RIT::ShankRIT) = RIT.Tw
transitiontemperature(RIT::ShankRIT) = RIT.Tc

"thermalpower(G, n, Tbath, T)
Return the thermal flow from a TES at temperature `T` and the bath at temperature `p.Tbath`.
Do the calculation carefully to avoid type instability from raising `T^n`."
function thermalpower(G, n, Tbath, T)
    g = Float64(G/u"W/K")
    t = Float64(T/u"K")
    tb = Float64(Tbath/u"K")
    k = g/n/t^(n-1) # (pW/K^n)
    power = k*(t^n-tb^n)
    return power*1u"W"
end
thermalpower(p::TESParams, T) = thermaflow(p.G, p.N, p.Tbath, T)





"Constructor that fixes `Tw` and `A` for a ShankRIT to have the given `alpha` and `beta`
parameters when biased at resistance `R0`."
function ShankRIT(alpha, beta, n, Tc, Tbath, G, R0, Rn)
    T0 = Tc /(1 + 3*beta/(2*alpha) - 2*(Rn-R0)/(Rn*alpha)*atanh(2*Float64(R0/Rn)-1))
    I0 = sqrt(thermalpower(G, n, Tbath, Tc)/R0)
    Tw = T0*(Rn-R0)/(Rn*log(2)*alpha)
    A = I0*(2*alpha/(3*T0*beta))^(3/2)
    ShankRIT(Tc,Rn,Tw,A)
end


type BiasedTES{T}
    p::TESParams{T}
    I0::typeof(1.0u"mA") # intial current for diff eq, aka current through TES
    T0::typeof(1.0u"mK") # initial temperature for diff equations, aka temperature of TES
    V ::typeof(1.0u"mV")  # Thévinen equivalent voltage V = I0*(p.Rl+p.R0), R0=quiescent resistance
                # also equal to Ibias*Rshunt. Careful! This V is a constant and is NOT
                # the voltage drop across the TES.
end


R(I,T, p::TESParams) = R(I,T, p.RIT)
"thermal TES equation
C*dT/dt = I^2*R - k*(T^n-Tbath^n)"
function dT(I, T, G, n, Tbath, C, R)
    # if you run into domain errors, try uncommenting the following line
    # but it probably represents a mistake in your timesteps that should be fixed.
    #T=max(T,0.0) # avoid domain errors when raising T to a power
    Q=I^2*R-thermalpower(G,n,Tbath,T)
    Q/C
end
"electrical TES equation
L*di/dt = (IBias-I)*Rs+I*Rs-I*Rtes"
function dI(I, T, V, Rl, L, R)
    (V-I*(Rl+R))/L
end
"Calling a BiasedTES gives the dI and dT terms for integration in an in place manner."
function (bt::BiasedTES)(t, u::Vector{Float64}, du::Vector{Float64})
    @show t,u,du
    T,I = max(0.0,u[1])*1u"K",u[2]*1u"A"
    p = bt.p
    r = R(I,T,p)
    du[1] = dT(I, T, p.G, p.n, p.Tbath, p.C, r)/u"K/s"
    du[2] = dI(I,T, bt.V, p.Rl, p.L, r)/u"A/s"
    @show du
    du
end

function pulses(nsample::Int, dt::Float64, bt::BiasedTES, Es::Vector, arrivaltimes::Vector; dtsolver=1e-9, method=DifferentialEquations.Tsit5(), abstol=1e-9, reltol=1e-9)
  u0 = Float64[bt.T0/u"K", bt.I0/u"A"]
  saveat = range(0,dt, nsample)
  prob = ODEProblem(bt, u0, (0.0, last(saveat)))
  Esdict = Dict([(at,E) for (at,E) in zip(arrivaltimes,Es)])
  # this defines a callback that is evaluated when t equals a value in arrival times
  # when evaluated it discontinuously changes the temperature (u[1])
  # the last (true,true) argument has to do with which points are saved
  function cbfun(integrator)
    integrator.u[1]+=Esdict[integrator.t]/bt.p.C/u"K"
    # modify the integrator timestep back to dtsolver, to take small steps on the rising edge of the pulse
    integrator.dtpropose=dtsolver # in future use modify_proposed_dt!, see http://docs.juliadiffeq.org/latest/basics/integrator.html#Stepping-Controls-1
  end
  cb = DiscreteCallback((t,u,integrator)->t in arrivaltimes, cbfun; save_positions=(false,false))
  # tstops is used to make sure the integrator checks each time in arrivaltimes
  sol = solve(prob,method,dt=dtsolver,abstol=abstol,reltol=reltol, saveat=saveat, save_everystep=false, dense=false,callback=cb, tstops=arrivaltimes)

  T = sol[1,:]*1u"K"
  I = sol[2,:]*1u"A"
  Rout = [R(I[i],T[i],bt.p) for i=1:length(T)]
  TESRecord(T,I, Rout,dt*1u"s")
end
type TESRecord
    T::Vector{typeof(1.0u"mK")}  # temperature (K)
    I::Vector{typeof(1.0u"mA")}  # current (A)
    R::Vector{typeof(1.0u"mΩ")}  # TES resistance (Ohm)
    dt::typeof(1.0u"s")  # seconds between samples (seconds)
end
# times(r::TESRecord) = range(0,r.dt,length(r.I))
Base.length(r::TESRecord) = length(r.I)
temperatures(r::TESRecord) = R.T
currents(r::TESRecord) = R.I
resistances(r::TESRecord) = R.R

end # module
