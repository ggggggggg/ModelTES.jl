using Revise
using ModelTES, Unitful, Test

# matches Kelsey's low-E pix
Tc=94.0u"mK"; Tbath = 75.0u"mK"
R0 = 1.64u"mΩ"; Rl = 0.35u"mΩ"; Rn = 8.2u"mΩ"; Rpara=0.0u"mΩ"
n=3.15; G = 27.34u"pW/K"; C = 0.1u"pJ/K"
alpha = 175.0; beta = 2.28; L=50.0u"nH"
rit = ShankRIT_from_αβ(alpha, beta, n, Tc, Tbath, G, R0, Rn);
tes_param = TESParams(n,Tbath,G,C,L,Rl,Rpara,rit)
I0=1.506454430784236e-2u"mA"
T0= 91.6286322033873u"mK"
V0=2.9978516417800104e-5u"mV"
bt0 = BiasedTES(tes_param, I0,T0,V0)

Teps = T0*1e-5
@test ModelTES.thermalpower(G, n, T0, T0+Teps) ≈ G*Teps rtol=1e-2

r = 1.4481096647342295u"mΩ"
@test rit(I0, T0) == r
@test ModelTES.dT(bt0.I0, bt0.T0, bt0.p.G, bt0.p.n, bt0.p.Tbath, bt0.p.C, r) ≈ -0.43425901199703454u"K/s"
@test ModelTES.dI(bt0.I0, bt0.T0, bt0.Vt, bt0.p.Rl, bt0.p.L, r) ≈ 0.057816274061034637u"A/s"
# bt(0,u,du)

bt1 = BiasedTES_from_R0(tes_param, R0)
r1 = ModelTES.rit(bt1)(bt1.I0, bt1.T0)
@test r1 ≈ R0
@test ModelTES.dT(bt1.I0, bt1.T0, bt1.p.G, bt1.p.n, bt1.p.Tbath, bt1.p.C, r1)|>u"K/s" ≈ 0u"K/s" atol = 1e-12u"K/s"
du = Float64[0, 0]
bt1(du, unitless.([bt1.T0/1u"K", bt1.I0/1u"A"]), nothing, nothing)
@test du ≈ [0,0] atol=1e-12


bt2 = BiasedTES(bt1.p, bt1.I0, bt1.T0, bt1.Vt)
bt2(du, unitless.([bt2.T0/1u"K", bt2.I0/1u"A"]), nothing, nothing)
out = pulses(1000, 10u"μs", bt2, [1000u"eV"], [100u"μs"], dtsolver=1u"μs")
plot(Float64.(unitless.(ModelTES.times(out)/1u"μs")), ModelTES.unitless.(out.T./out.T[1]), label="T/T0")
plot(Float64.(unitless.(ModelTES.times(out)/1u"μs")), ModelTES.unitless.(out.I./out.I[1]), label="I/I0")
plot(Float64.(unitless.(ModelTES.times(out)/1u"μs")), ModelTES.unitless.(out.R./out.R[1]), label="R/R0")
xlabel("time (μs)")
ylabel("T (mK)")
legend()
ylim(0.5,1.5)
;
# using ModelTES, ARMA
# using Test

# biased_tess = [ModelTES.pholmes(), ModelTES.lowEpix(), ModelTES.highEpix(), ModelTES.LCLSII()]

# for bt in biased_tess
#   Vs_in = bt.V*collect(0:0.1:10)
#   Is, Ts, Rs, Vs_out = ModelTES.iv_curve(bt.p, Vs_in)
#   derivs = zeros(Float64, length(Vs_in), 2)
#   for i in eachindex(Vs_in)
#       derivs[i,:] = ModelTES.dT_and_dI_iv_point(bt.p, Is[i], Ts[i], Rs[i], Vs_in[i])
#   end
#   @test maximum(abs.(derivs))<1e-10

#   #compare rk8 and DifferentialEquations intergrators
#   out = rk8(12000,1e-7, bt, 1000, 2000);
#   out_pulse = pulse(12000,1e-7, bt, 1000, 2000);
#   function worst_relative_error(a,b)
#       @assert(all(times(a).==times(b)))
#       eI=maximum(abs.(2*(a.I-b.I)./(a.I.+b.I)))
#       eT=maximum(abs.(2*(a.T-b.T)./(a.T.+b.T)))
#       eR=maximum(abs.(2*(a.R-b.R)./(a.R.+b.R)))
#       max(eI,eT,eR)
#   end
#   @test worst_relative_error(out,out_pulse)<1e-5
#   out_temp = deepcopy(out)
#   out_temp.I[1]*=1.11
#   @test worst_relative_error(out_temp,out_pulse)>1e-1
#   out_temp = deepcopy(out)
#   out_temp.T[1]*=1.11
#   @test worst_relative_error(out_temp,out_pulse)>1e-1
#   out_temp = deepcopy(out)
#   out_temp.T[1]*=1.11
#   @test worst_relative_error(out_temp,out_pulse)>1e-1

#   # compare to a pulse output with bigger timesteps, adapative solving should make this work
#   out_for_resample = rk8(12000,1e-7, bt, 1000,0);
#   out_pulse_ts = pulse(1200,1e-6, bt, 1000,0);
#   out_ts = ModelTES.TESRecord(out_for_resample.T[1:10:end], out_for_resample.I[1:10:end], out_for_resample.R[1:10:end], 1e-6)
#   @test worst_relative_error(out_ts,out_pulse_ts)<1e-5
#   out_pulse_ts2 = pulse(120,1e-5, bt, 1000,0);
#   out_ts2 = ModelTES.TESRecord(out_for_resample.T[1:100:end], out_for_resample.I[1:100:end], out_for_resample.R[1:100:end], 1e-5)
#   @test worst_relative_error(out_ts2,out_pulse_ts2)<1e-5


#   # Integrate a pulse with 12000 samples, 1e-7 second spacing, 1000 eV energy, 2000 presamples from the higher biased version of the same tes
#   out2 = pulse(12000,1e-7, bt, 1000, 2000);


#   # many pulses in one trace
#   outmany = ModelTES.pulses(12000,1e-7, bt, [1000,1000,2000,3000,1000,500,2000], collect(1:7)*2e-4);
#   @test length(times(outmany))==length(outmany.I)
#   #make the pulses arrive halfway between time points
#   outmany2 = ModelTES.pulses(12000,1e-7, bt, [1000,1000,2000,3000,1000,500,2000], 0.5e-7 .+ collect(1:7)*2e-4);
#   @test times(outmany)==times(outmany2)
#   @test length(outmany)==length(outmany2)==12000

#   # compare the difference between when the pulses arrive half way between time points, and 1 time point apart
#   # the integrated difference should be about a factor of two apart
#   a=sum(abs.(outmany.I[2:end-1]-outmany2.I[2:end-1])) # pulses off by half a sample
#   b=sum(abs.(outmany.I[2:end]-outmany.I[1:end-1])) # pulses off by one sample
#   @test isapprox(a,b/2,rtol=1e-2,atol=1e-5)

#   # Compute the noise spectrum etc.
#   nmodel = NoiseModel(bt, 1e-6)
#   freq = range(0, stop=5e5, length=50)
#   psd = model_psd(nmodel, freq)
#   covar = model_covariance(nmodel, 20)

# end

# include("ihtes.jl")
