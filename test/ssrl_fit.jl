using Revise
using ModelTES, Unitful, Test, PyPlot, DelimitedFiles, PrettyTables
Revise.includet("lmfit.jl")

pulse_types = ["FeLa", "CaKa", "PdLa", "SiKa"]
energies_ev = [704.8, 3691.719, 2838.638, 1739.985]
paths = ["/Users/oneilg/Documents/pulse fitting/ssrl_pixel_data/20200212_B3all_0p95V_FeLa.dat",
"/Users/oneilg/Documents/pulse fitting/ssrl_pixel_data/20200214_B3all_0p95V_CaKa.dat",
"/Users/oneilg/Documents/pulse fitting/ssrl_pixel_data/20200214_B3all_0p95V_PdLa.dat",
"/Users/oneilg/Documents/pulse fitting/ssrl_pixel_data/20200214_B3all_0p95V_SiKa.dat",
]
tss = Vector{Float64}[]
pulses = Vector{Float64}[]
for (pt, path) in zip(pulse_types, paths)
    data = readdlm(path)
    ts = float.(data[2:end,1])
    pulse = float.(data[2:end,2])
    push!(tss, ts)
    push!(pulses, pulse)
end

figure()
for (i,(ts, pulse)) in enumerate(zip(tss, pulses))
    plot(ts, pulse, label="$(energies_ev[i])")
end
legend()
xlabel("times (s)")
ylabel("pulse (arbs)")

# for fitting we want one concatenated pulse
ydata = vcat(pulses...)
dt_s = tss[1][2]-tss[1][1]
dt = dt_s*1u"s"
xdata = range(0, step=dt_s, length=length(ydata))

figure()
plot(xdata, ydata)
xlabel("times (s)")
ylabel("pulse (arbs)")

first_pulse_arrival_index = 514.5
first_pulse_arrival_s = first_pulse_arrival_index*dt_s
pulse_arrivals = collect(range(first_pulse_arrival_s, step=dt_s*length(pulses[1]), length=4).*1u"s")
pulse_energies = energies_ev.*1u"eV"

# define our starting guesses
Rshunt = 234.4e-3u"mΩ";
Rbias = 9280u"Ω"; Vbias = 0.95u"V"; expected_Vt = Vbias*Rshunt/Rbias |> u"μV";
Tc=60.7u"mK"; Tbath = 52.0u"mK"
R0 = 1.0u"mΩ"; Rpara=1e-3u"mΩ"; Rl = Rshunt+Rpara; Rn = 9.29u"mΩ"; 
n=4.03; G = 33.3u"pW/K"; C = 0.3u"pJ/K"; L=50.0u"nH"
target_alpha = 400
target_beta = 1
rit = ShankRIT_from_αβ(target_alpha, target_beta, n, Tc, Tbath, G, R0, Rn)
tes_param = TESParams(n,Tbath,G,C,L,Rl,Rpara,rit)
bt0 = BiasedTES_from_R0(tes_param, R0)
Ic0 = 1u"mA"
rit2 = ModelTES.TwoFluidRIT_from_α(target_alpha, R0, bt0.T0, bt0.I0, Ic0, Rn, Tc)
tes_param2 = TESParams(n,Tbath,G,C,L,Rl,Rpara,rit2)
bt2 = BiasedTES_from_R0(tes_param2, R0)

out = ModelTES.pulses(length(ydata), dt, bt2, pulse_energies, pulse_arrivals, dtsolver=1u"μs")
arbs_per_current = -300.0u"μA^-1" 
arbs = unitless.((out.I.-bt2.I0)*arbs_per_current)

figure()
plot(xdata, arbs, label="first guess")
plot(xdata, ydata, label="data")
xlabel("time (s)")
ylabel("pulse in arbs")

# define fitting helper functions
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
function make_2fluid_biased_tes(n, Tbath, G, C, L, Rl, Rp, ci, Rn, Ic0, Tc, I0, T0, Vt, _arbs_per_current)
  rit = ModelTES.TwoFluidRIT(ci, Rn*1u"mΩ", Ic0*1u"mA", Tc*1u"mK")
  p = ModelTES.TESParams(n, Tbath*1u"mK", G*1u"pW/K",
      C*1u"pJ/K", L*1u"nH", Rl*1u"mΩ", Rp*1u"mΩ",  rit)
  ModelTES.BiasedTES(p, I0*1u"mA", T0*1u"mK", Vt*1u"mV")
end  
try # allows redefinition of fitfunc more easily
    for m in methods(fitfunc) 
        Base.delete_method(m) 
    end
catch
end
function fitfunc(x, n, Tbath, G, C, L, Rl, Rp, ci, Rn, Ic0, Tc, I0, T0, Vt, _arbs_per_current)
  bt = make_2fluid_biased_tes(n, Tbath, G, C, L, Rl, Rp, ci, Rn, Ic0, Tc, I0, T0, Vt)
  out = ModelTES.pulses(length(ydata), dt, bt, pulse_energies, pulse_arrivals, dtsolver=1u"μs")
  arbs_per_current = _arbs_per_current*1u"μA^-1"
  arbs = unitless.((out.I.-bt.I0)*arbs_per_current)
  return arbs
end
model, params = model_and_params(fitfunc);
copy_fields!(params, bt2);
params.n(vary=false, min=3, max=5)
params.G(vary=true, min=1, max=1000)
params.Rl(vary=false)
params.Tbath(vary=true, min=30, max=100)
params.L(vary=true, min=25, max=1e4)
params.Rn(vary=false)
params.Rp(vary=false)
params.Vt(vary=false)
params.ci(val=0.7, vary=true)
params._arbs_per_current(val=-351, vary=false)

weights = float.(xdata .< 4000);
result = fit(model, params, x=xdata, y=ydata, weights=weights);

figure(figsize=(14,6))
plot(xdata, ydata, label="data")
plot(xdata, result(x=xdata, params=result.params), label="fit")
plot(xdata, result(x=xdata, params=result.init_params), label="guess")
xlabel("time (s)")
ylabel("signal (arbs)")
legend()
plt.tight_layout()
bt_fit = make_2fluid_biased_tes(values(result.params)...)
ModelTES.alpha_beta(bt_fit)
;

function ptable(result)
pretty_table((name=[p.name for p in result.params],
unit=["", "mK", "pW/K", "pJ/K", "nH", "mΩ", "mΩ", "", "mΩ", "mA", "mK", "mA", "mK", "V", "μA^-1"],
vary=[p.vary for p in result.params],
guess=values(result.init_params),
value=values(result.params),
fit_uncertainty=[p.unc==nothing ? NaN : p.unc for p in result.params],
), 
formatters = ft_printf("%5.3g", [3,4,5,6]))
end
ptable(result)
