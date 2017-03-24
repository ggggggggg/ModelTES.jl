# Unvalidated code for attempting a stochastic integration of the differential
# equations for the TES. Galen's project.




function stochasticdy!(dy,y,dt,bt::BiasedTES)
   p=bt.p
   T, I = y # tes tempertature, tes current
   r = R(I,T,p) #res resistance
   # generate random variables for use in stocastic integration
   bath_thermal_noise, tes_thermal_noise, rl_johnson_noise, tes_johnson_noise = randn(4)

   johnsontes = tes_johnson_noise*2*sqrt(kb*T*r/dt)
   johnsonrl = rl_johnson_noise*2*sqrt(kb*p.Tbath*p.Rl/dt)

   # k*T^n*dt/kb*T is the average number of phonons going from the TES to the bath
   # assuming all phonons of energy kb*t
   # sqrt(k*T^n*dt/kb*T) is the std deviation of the number of phonons
   # the same reasoning applies in the other direction
   dqdt = -p.k*(T^p.n-p.Tbath^p.n)+I^2*r +
     tes_thermal_noise*sqrt(kb*p.k/dt)*T^(p.n/2 + 0.5) + # variation in number of phonon TES->bath
     bath_thermal_noise*sqrt(kb*p.k/dt)*p.Tbath^(p.n/2+0.5) + # variation in number of phonon bath->TES
     -I*johnsontes # work on the noise voltage source for the TES add/remove energy from the TES

   # the resistors are modeled as noiseless resistors in parallel with
   # voltage sources having white noise
   # the load resistor is assumed to be at Tbath
   ldidt = bt.V-I*(p.Rl+r) + johnsontes + johnsonrl

   if isnan(dqdt)
      @show dy, y, dt, bt
      error()
   end


   dy[1] = dt*dqdt/p.C
   dy[2] = dt*ldidt/p.L
end


"stochastic(nsample::Int, dt::Float64, bt::BiasedTES, E::Number, npresample::Int=0, nrandomizingsample::Int=6000)
Integrate a pulse record with with stocastic noise `E` eV with `nsample` total samples, `npresample` presamples, and use
`nrandomizingsample` samples thrown away to ensure randomized starting `T0` and `I0`. Use `E=0` for noise records."
function stochastic(nsample::Int, dt::Float64, bt, E::Number, npresample::Int=0, nrandomizingsample::Int=6000)
   p = bt.p
   T = Array(Float64, nsample)
   I = Array(Float64, nsample)
   y = [bt.T0, bt.I0]
   dy = Array(Float64,2)
   # manually implement euler's method for integration
   for i = 1:nrandomizingsample # integrate for a while to get well randomized T0 and I0
       stochasticdy!(dy,y,dt,bt)
       y+=dy
   end
   for i = 1:nsample
      if i == npresample+1
         y[1]+=E*J_per_eV/p.C
      end # increase temperature to start pulse
      stochasticdy!(dy,y,dt,bt)
      y+=dy
      T[i] = y[1]
      I[i] = y[2]
   end
   TESRecord(T,I,R(I,T,bt.p),dt)
end
