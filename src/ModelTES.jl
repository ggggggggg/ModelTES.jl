module ModelTES
export
    BiasedTES,
    BiasedTES_from_R0,
    TESParams,
    ShankRIT,
    ShankRIT_from_αβ,
    pulses,
    unitless

using Roots, ForwardDiff, DifferentialEquations, Unitful
"""assert `x` is unitless and return a Number, eg Float64 or Rational{Int}"""
unitless(x) = uconvert(Unitful.NoUnits, x)


abstract type AbstractRIT end

# following Irwin-Hilton figure 3
struct TESParams{RITType<:AbstractRIT}
    n       ::Float64   # thermal conductance exponent (unitless)
    Tbath   ::typeof(1.0u"mK")  # bath temperature

    G       ::typeof(1.0u"pW/K")  # thermal conductivity G = n*K*(Tc^n-1)
    C       ::typeof(1.0u"pJ/K")   # heat capacity of TES

    L       ::typeof(1.0u"nH") # inductance of SQUID (H)
    Rl      ::typeof(1.0u"mΩ")  # Thevenin-equivalent resistance Rload = (Rshunt+Rparasitic)(ohms)
                        # note Ibias = V/Rshunt
    Rp      ::typeof(1.0u"mΩ")    # Rparastic; Rshunt = Rl-Rparasitic (ohms)

    RIT     ::RITType   # RIT surface

    _K       ::Float64  # W/K^n, dont want type instability with n, so we use a Float64
                        # but we do want to pre-compute this so we can avoid some ^n calculations
                        # redundant with G, so make sure they are consistent!
end
rit(p::TESParams) = p.RIT
"""    _K_from_G(G, Tc, n)
return `K` as a Float64 in units W/K^n, used to avoid types depending on n
"""
function _K_from_G(G, Tc, n)
    Tc = unitless(Tc/u"K")
    _K = unitless(G/u"W/K")/(n*Tc^(n-1)) # W/K^n
end    
"""    TESParams(n, Tbath, G, C, L, Rl, Rp, RIT)
Define TESParams, does unit conversion, and calculates _K from G"""
function TESParams(n, Tbath, G, C, L, Rl, Rp, RIT)
    Tc = transitiontemperature(RIT)
    _K = _K_from_G(G, Tc, n)
    TESParams(n, Tbath |> u"mK", G |> u"pW/K", C |> u"pJ/K", L |> u"nH", Rl |> u"mΩ", Rp |> u"mΩ", RIT, _K)
end
struct BiasedTES{RITType}
    p::TESParams{RITType}
    I0::typeof(1.0u"mA") # intial current for diff eq, aka current through TES
    T0::typeof(1.0u"mK") # initial temperature for diff equations, aka temperature of TES
    Vt ::typeof(1.0u"mV")  # Thévinen equivalent voltage Vt = I0*(p.Rl+p.R0), R0=quiescent resistance
                # also equal to Ibias*Rshunt. 
                # Careful!  Vt is a constant and is NOT the voltage drop across the TES.
end
function BiasedTES(p, I0, T0, Vt)
    BiasedTES(p, uconvert(u"mA", I0), uconvert(u"mK", T0), uconvert(u"mV", Vt))
end
rit(bt::BiasedTES) = rit(bt.p)

function Base.show(io::IO, p::T) where T <: Union{TESParams, AbstractRIT, BiasedTES}
    print(io, T, "(")
    for (i,fname) in enumerate(fieldnames(T))
        print(io, "$fname=", getfield(p, fname), i<fieldcount(T) ? ", " : "")
    end
    print(")")
end


"Find the initial conditions (I0, T0, V) that cause `p` to have resistance `targetR`."
function initialconditions(p::TESParams, targetR)
    Tc = transitiontemperature(rit(p))
    Tw = transitionwidth(rit(p))
    # we want to find a combination of I and T such that
    # 1. R(I,T) = targetR
    # 2. the steady state temperature of the TES is T
    # The steady state temperature occurs when Pjoule = Pthermal
    # where Pjoule = I^2*R
    # and Pthermal = thermalpower(p, T) = k*(T^n-Tbath^n)
    "For a given T0, the difference R-targetR. Use to solve numerically for T0."
    function getR0error(T0)
        I0 = sqrt(thermalpower(p, T0)/targetR)
        rit(p)(I0,T0)-targetR
    end
    T00 = fzero(getR0error, Tc-10*Tw, Tc+10*Tw)
    I00 = sqrt(thermalpower(p, T00)/targetR) |> u"nA"
    R00 = rit(p)(I00,T00)
    V00 = I00*(p.Rl+R00)
   # now evolve these conditions through integration to really lock them in.
   # shouldn't hard code step size here
    if false # turn off the diff eq part
        out = pulses(10, 10u"ms", BiasedTES(p, I00, T00, V00), [0u"eV"], [0u"s"])
        T0 = out.T[end]
        I0 = out.I[end]
        R0 = rit(p)(I0, T0)
        V = I0*(p.Rl+R0)
        return I0, T0, V
    else
        return I00, T00, V00
    end
end


"Created a biased tes with quiescent state resistance R0"
function BiasedTES_from_R0(p::TESParams{RITType}, R0) where RITType
   @assert 0u"Ω" < R0 < normal_resistance(rit(p))
   I0, T0, V = initialconditions(p, R0)
   BiasedTES(p, I0, T0, V)
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
    Is=Vector{Float64}(undef, length(Vs))
    Ts=Vector{Float64}(undef, length(Vs))
    Rs=Vector{Float64}(undef, length(Vs))
    Vs_out=Vector{Float64}(undef, length(Vs))
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



"Calculate `R0` the quiescent resistance of `bt`."
getR0(bt::BiasedTES) = rit(bt)(bt.I0, bt.T0)
getG0(bt::BiasedTES) = bt.p.G




"Calculate paramaters in Irwin-Hilton table 1."
function getlinearparams(bt::BiasedTES)
   f(x) = rit(bt)(x[1], x[2])
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
   bt.I0, bt.T0, bt.Vt, p.Rl, p.Tbath, p.Tbath, p.L, R0, G0, p.C, alpha, beta, loopgain,
        tauthermal, taucc, taueff, tauelectrical, tauplus, tauminus, lcritplus, lcritminus, lcritical
end

"Paramters from Irwin-Hilton table one for modeling a linear TES. Defined in Table 1 of Irwin-Hilton chapter."
struct IrwinHiltonTES
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

struct ShankRIT <: AbstractRIT
    Tc::typeof(1.0u"mK")
    Rn::typeof(1.0u"mΩ")
    Tw::typeof(1.0u"mK")  # transition width
    A::typeof(1.0u"A/K^(3/2)")  # current dependence for R(T,I) (A/K^(3/2))
end
"R(I, T, RIT::ShankRIT)
Return resistance as a function of current `I` and temperature `T` using the model `RIT`.
TES resistance model (Shank et al. 2014)"
function (RIT::ShankRIT)(I, T)
    z = cbrt(I/RIT.A)
    dTc = z*z
    # unitful and the ^ operator dont play super nice
    # so use cbrt and squaring by multiplication so do a 2/3 power
    # dTc = (I/A)^(3/2)
    x = unitless((T-RIT.Tc+dTc)/(2*log(2)*RIT.Tw))
    RIT.Rn/2*(1+tanh(x))
end
transitionwidth(RIT::ShankRIT) = RIT.Tw
transitiontemperature(RIT::ShankRIT) = RIT.Tc
normal_resistance(RIT::ShankRIT) = RIT.Rn

"Constructor that fixes `Tw` and `A` for a ShankRIT to have the given `alpha` and `beta`
parameters when biased at resistance `R0`."
function ShankRIT_from_αβ(alpha, beta, n, Tc, Tbath, G, R0, Rn)
    T0 = Tc /(1 + 3*beta/(2*alpha) - 2*(Rn-R0)/(Rn*alpha)*atanh(2*Float64(R0/Rn)-1))
    _K = _K_from_G(G, Tc, n)
    I0 = sqrt(thermalpower(_K, n, Tbath, T0)/R0)
    Tw = T0*(Rn-R0)/(Rn*log(2)*alpha)
    A = I0*(2*alpha/(3*T0*beta))^(3/2)
    ShankRIT(Tc,Rn,Tw,A)
end

"""
We follow Bennet and Ullom 2015 equations 33 and 36. 
We notice that `cr` and `Rn` only ever appear as `cr*Rn`, so we define cr=1 and simply
do not include cr in the code.
Morgan 2017 provides some values for the parameters of the model for real devices.
    eg. ci=0.85, Ic0 ~200*I0
"""
struct TwoFluidRIT <: AbstractRIT
    ci::Float64
    Rn::typeof(1.0u"mΩ")
    Ic0::typeof(1.0u"mA")
    Tc::typeof(1.0u"mK")
end
function (tf::TwoFluidRIT)(I,T)
    # following Bennet and Ullom 2015
    t = clamp(unitless(T/tf.Tc), 0.0, 1.0) # avoid domain error in next line, also avoid exponentiation unitful quantities
    Ic = tf.Ic0*(1-t)^(3//2) # eq. 33
    if I < tf.ci*Ic # per Doug, there is only resistance if I>ci*Ic
        return 0u"mΩ" # ensure both branches return the exact same units
    else
        return tf.Rn*(1-tf.ci*Ic/I) |> u"mΩ" # eq. 36
    end
end
transitiontemperature(tf::TwoFluidRIT) = tf.Tc
transitionwidth(tf::TwoFluidRIT) = tf.Tc/100 # hack to move forward, the correct answer is I0 dependent, maybe define this based on 1 nA?
normal_resistance(tf::TwoFluidRIT) = tf.Rn

function TwoFluidRIT_from_α(alpha, R0, T0, I0, Ic0, Rn, Tc)
    b = (3/2)*(Rn/R0)*(Ic0/I0)*(T0/Tc)*sqrt(1-T0/Tc)
    ci = alpha/b
    return TwoFluidRIT(ci, Rn, Ic0, Tc)
end

"""    thermalpower(K, n, Tbath, T)
`K` is a Float64, representing a number with units `W/K^n`. `Tbath` and `T` should have temperature units.
Return the thermal flow from a TES at temperature `T` and the bath at temperature `p.Tbath`.
Do the calculation carefully to avoid type instability from raising `T^n`."""
function thermalpower(K, n, Tbath, T)
    t = unitless(T/u"K")
    tb = unitless(Tbath/u"K")
    power = K*(t^n-tb^n)
    return power*1u"W"
end

"`Z(tes::IrwinHiltonTES, f)`

Returns the impedance of the `tes` at frequency `f`.
Implements equation 42 of Irwin-Hilton chapter."
function Z(tes::IrwinHiltonTES, f)
  ω=2π*f
  tes.R0*(1+tes.beta) .+ tes.R0*tes.loopgain*(2+tes.beta)./((1-tes.loopgain)*(1 .+im*ω*tes.taucc))
end
thermalpower(p::TESParams, T) = thermalpower(p._K, p.n, p.Tbath, T)

"`Zcircuit(tes::IrwinHiltonTES, f)`

Returns impedance of complete circuit of `tes` at frequency `f`."
function Zcircuit(tes::IrwinHiltonTES, f)
  ω=2π*f
  tes.R0 .+ im*ω*tes.L .+ Z(tes,ω)
end



"thermal TES equation
C*dT/dt = I^2*R - k*(T^n-Tbath^n)"
function dT(I, T, K, n, Tbath, C, R)
    # if you run into domain errors, try uncommenting the following line
    # but it probably represents a mistake in your timesteps that should be fixed.
    # T=max(T,0.0u"K") # avoid domain errors when raising T to a power
    Q=I^2*R-thermalpower(K,n,Tbath,T)
    Q/C
end
"electrical TES equation
L*di/dt = (IBias-I)*Rs+I*Rs-I*Rtes
where Vt = (IBias-I)*Rs+I*Rs = IBias*Rs"
function dI(I, T, Vt, Rl, L, R)
    (Vt-I*(Rl+R))/L
end

"""Make BiasedTES callable as for differntial equation solver"""
function (bt::BiasedTES)(du, u, p_, t) # use DifferentialEquation 4.0+ API
    T,I = max(0.0,u[1])*1u"K",u[2]*1u"A"
    # @show round(t*1e6), T, I
    p = bt.p
    r = rit(bt)(I,T)
    du[1] = unitless(dT(I, T, p._K, p.n, p.Tbath, p.C, r)/u"K/s")
    du[2] = unitless(dI(I,T, bt.Vt, p.Rl, p.L, r)/u"A/s")
    du
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
    T = Array{Float64}(undef, nsample)
    I = Array{Float64}(undef, nsample)
    T[1:npresamples] .= bt.T0
    I[1:npresamples] .= bt.I0 # set T0, I0 for presamples
    y = [bt.T0+E*J_per_eV/p.C, bt.I0]; ys = similar(y); work = Array{Float64}(undef, 14)
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
    Rout = [rit(bt)(I[i],T[i]) for i=1:length(T)]

    TESRecord(T, I, Rout, dt)
end

"""    pulses(nsample::Int, dt, bt::BiasedTES, Es::Vector, arrivaltimes::Vector; 
    dtsolver=1e-9, method=DifferentialEquations.Tsit5(), abstol=1, reltol=1e-7)

return a record with photons arriving at `arrivaltimes` with energies `Es`"""
function pulses(nsample::Int, dt, bt::BiasedTES, Es, arrivaltimes; 
        dtsolver=100u"ns", method=DifferentialEquations.Rodas4P(), abstol=1, reltol=1e-7)
    # convert to unitless numbers
    u0 = Float64.(unitless.([bt.T0/1u"K", bt.I0/1u"A"])) # these have different units! we want u to be a Vector of a single type
    dt = Float64(unitless(dt/1u"s"))
    arrivaltimes = Float64.(unitless.(arrivaltimes./1u"s"))
    dtsolver = Float64(unitless(dtsolver/1u"s"))
    saveat = range(0, step=dt, length=nsample)
    prob = ODEProblem(bt, u0, (0.0, last(saveat)))
    @assert length(Es) == length(arrivaltimes)
    i = 1 # used in callback
    n = length(Es)
    # this defines a callback that is evaluated when t equals a value in arrival times
    # when evaluated it discontinuously changes the temperature (u[1])
    # the last (true,true) argument has to do with which points are saved
    function cbfun(integrator)
        integrator.u[1]+=unitless(Es[i]/bt.p.C/1u"K")
        i+=1
        # modify the integrator timestep back to dtsolver, to take small steps on the rising edge of the pulse
        set_proposed_dt!(integrator, dtsolver)
    end
    cb = DiscreteCallback((u,t,integrator)->(t in arrivaltimes), cbfun; save_positions=(false, false))
    # as of Jul 2020, the DiffEq docs say discontinuously changing u requires using save before and after
    # for now im ignoring this and hoping `set_proposed_dt!` is sufficient
    # tstops is used to make sure the integrator checks each time in arrivaltimes
    sol = solve(prob, method=method, dt=dtsolver, abstol=abstol, reltol=reltol, saveat=saveat, 
        save_everystep=false, dense=false, callback=cb, tstops=arrivaltimes, adaptive=false)

    v = sol(saveat)
    T = v[1,:]*1u"K"
    I = v[2,:]*1u"A"
    Rout = [rit(bt)(I[i],T[i]) for i=1:length(T)]
    TESRecord(T,I, Rout,dt*1u"s")
end

struct TESRecord
    T::Vector{typeof(1.0u"mK")}  # temperature (K)
    I::Vector{typeof(1.0u"mA")}  # current (A)
    R::Vector{typeof(1.0u"mΩ")}  # TES resistance (Ohm)
    dt::typeof(1.0u"s")  # time between samples (seconds)
end
times(r::TESRecord) = range(0u"s", step=r.dt, length=length(r.I))
Base.length(r::TESRecord) = length(r.I)
temperatures(r::TESRecord) = R.T
currents(r::TESRecord) = R.I
resistances(r::TESRecord) = R.R
function Base.show(io::IO, r::T) where T <: TESRecord
    print(io, T, "(")
    print(io, "length=$(length(r))")
    if length(r) >= 1
        print(io, ", T[1]=", r.T[1], ", I[1]=", r.I[1], ", R[1]=", r.R[1])
    end
    print(io, ", dt=", r.dt,")")
end


end # module
