using Revise
using ModelTES, Unitful, Test, PyPlot
Revise.includet("lmfit.jl")
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

Teps = 1u"μK"
@test ModelTES.thermalpower(tes_param._K, n, Tc, Tc+Teps) ≈ G*Teps rtol=1e-2

r = 1.7467637973165189u"mΩ"
@test rit(I0, T0) == r
@test ModelTES.dT(bt0.I0, bt0.T0, bt0.p._K, bt0.p.n, bt0.p.Tbath, bt0.p.C, r) |> u"K/s" ≈ 0.44238556587270167u"K/s"
@test ModelTES.dI(bt0.I0, bt0.T0, bt0.Vt, bt0.p.Rl, bt0.p.L, r)|>u"A/s" ≈ -0.03216549419908801u"A/s"
# bt(0,u,du)

bt1 = BiasedTES_from_R0(tes_param, R0)
r1 = ModelTES.rit(bt1)(bt1.I0, bt1.T0)
@test r1 ≈ R0
@test ModelTES.dT(bt1.I0, bt1.T0, bt1.p._K, bt1.p.n, bt1.p.Tbath, bt1.p.C, r1)|>u"K/s" ≈ 0u"K/s" atol = 1e-12u"K/s"
du = Float64[0, 0]
bt1(du, unitless.([bt1.T0/1u"K", bt1.I0/1u"A"]), nothing, nothing)
@test du ≈ [0,0] atol=1e-12


Ic0 = 200 * bt1.I0
rit2 = ModelTES.TwoFluidRIT_from_α(alpha, R0, bt1.T0, bt1.I0, Ic0, Rn, Tc)
tes_param2 = TESParams(n,Tbath,G,C,L,Rl,Rpara,rit2)
bt2 = BiasedTES_from_R0(tes_param2, R0)

Ts = 90u"mK":.01u"mK":95u"mK"
figure()
R = rit2.(I0, Ts)
plot(unitless.(Ts./u"mK"), unitless.(R./u"mΩ"), label="I0 in title")
R = rit2.(1.3*I0, Ts)
plot(unitless.(Ts./u"mK"), unitless.(R./u"mΩ"), label="I0=I0_title*1.3")
Rshank = rit.(I0, Ts)
plot(unitless.(Ts./u"mK"), unitless.(Rshank./u"mΩ"), label="shank")
xlabel("T (mk)")
ylabel("R (mΩ)")
title("cr=1, ci=$(rit2.ci), I0=Ic0/200 (blue), Rn=8.2mΩ")
legend()

println("2fluid alpha_beta: ", ModelTES.alpha_beta(bt2))
println("shank alpha_beta: ", ModelTES.alpha_beta(bt1))


es = [500, 1000, 2000, 3000, 4000].*1u"eV"
tpulse = 4000u"μs" 
arrivals = range(0u"μs", step=tpulse, length=length(es))
dt = 5u"μs"
steps = round(Int, ModelTES.unitless((arrivals[end]+tpulse)/dt)) 

out1 = pulses(steps, dt, bt1, es, arrivals, dtsolver=1u"μs")
out2 = pulses(steps, dt, bt2, es, arrivals, dtsolver=1u"μs")

figure()
plot(Float64.(unitless.(ModelTES.times(out1)/1u"μs")), ustrip.(u"mA", out1.I), label="Shank")
plot(Float64.(unitless.(ModelTES.times(out2)/1u"μs")), ustrip.(u"mA", out2.I), label="2Fluid")
xlabel("time (μs)")
ylabel("I (mA)")
legend()
title(prod("$e " for e in es))
;close("all");




steps = 4000
# make data with ShankRIT
out_data = pulses(steps, dt, bt1, es, arrivals, dtsolver=1u"μs")
ydata = ustrip(out_data.I)
xdata = Float64.(unitless.(ModelTES.times(out_data)/1u"μs"))

function copy_fields!(params, x)
  for fname in fieldnames(typeof(x))
    if string(fname) in keys(params)
      params[string(fname)](val=ustrip(getproperty(x, fname)))
    else 
      copy_fields!(params, getproperty(x, fname))
    end
  end
  return params
end

# define a fitting model
function make_2fluid_biased_tes(n, Tbath, G, C, L, Rl, Rp, ci, Rn, Ic0, Tc, I0, T0, Vt)
  rit = ModelTES.TwoFluidRIT(ci, Rn*1u"mΩ", Ic0*1u"mA", Tc*1u"mK")
  p = ModelTES.TESParams(n, Tbath*1u"mK", G*1u"pW/K",
      C*1u"pJ/K", L*1u"nH", Rl*1u"mΩ", Rp*1u"mΩ",  rit)
  ModelTES.BiasedTES(p, I0*1u"mA", T0*1u"mK", Vt*1u"mV")
end  
function fitfunc(x, n, Tbath, G, C, L, Rl, Rp, ci, Rn, Ic0, Tc, I0, T0, Vt)
  bt = make_2fluid_biased_tes(n, Tbath, G, C, L, Rl, Rp, ci, Rn, Ic0, Tc, I0, T0, Vt)
  out = pulses(steps, dt, bt, es, arrivals, dtsolver=1u"μs")
  return ustrip(out.I)
end
model, params = model_and_params(fitfunc)
params.Tbath(vary=false)
params.L(vary=false)
params.Rn(vary=false)
params.Rp(vary=false)
copy_fields!(params, bt2)



result = fit(model, params, x=xdata, y=ydata)

ydata_init = result(params=result.init_params)
ydata_fit = result()
bt_fit = make_2fluid_biased_tes(values(params)...)


figure()
plot(xdata, ydata, label="fake data (from shank rit)", lw=2)
plot(xdata, ydata_init, label="initial guess (2fluid)", lw=2)
plot(xdata, ydata_fit, label="lsq fit (2fluid)", lw=2)
legend()
xlabel("time (μs)")
ylabel("current  (A)")

figure()
plot(xdata, ydata, label="fake data (from shank rit)", lw=2)
plot(xdata, ydata_fit, label="lsq fit (2fluid)", lw=2)
legend()
xlabel("time (μs)")
ylabel("current  (A)")

println("2fluid alpha_beta (init): ", ModelTES.alpha_beta(bt2))
println("shank alpha_beta: ", ModelTES.alpha_beta(bt1))
println("2fluid alpha_beta (fit): ", ModelTES.alpha_beta(bt_fit))