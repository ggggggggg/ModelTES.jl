using Revise
using ModelTES, Unitful, Test, PyPlot, DelimitedFiles, PrettyTables
using HDF5
Revise.includet("lmfit.jl")

rns = [10, 15, 20, 25, 30, 40]
h5 = h5open("pulses_vs_bias.hdf5", "r")
pulses_by_rn = Dict(rn => read(h5["$rn"]) for rn in rns) 
energies_kev = h5["energies"] |> read
times = h5["times"] |> read
raw_data_n_samples = size(pulses_by_rn[10], 1)
close(h5)

figure()
for (i, rn) in enumerate(rns)
    pulses = pulses_by_rn[rn]
    plot(times .+ (i - 1) * (times[end] - times[1]), pulses * 1e6)
end
legend()
xlabel("times (s)")
ylabel("pulse current (uA)")
;
# for fitting we want:
#     one concatenated pulse
#     some level of decimation
#     select some subset of pulses
#     a corresponding list of rn and energy
dt_s_no_decimation = times[2] - times[1]
decimation_factor = 10
n_samples_pulse = ceil(Int, raw_data_n_samples / decimation_factor)
dt_s = dt_s_no_decimation * decimation_factor
use_rns = [10, 25, 40]
# use_rns = [10]
n_rns = length(use_rns)
use_energy_inds = [1, 5, 10, 15]
# use_energy_inds = [5]
n_energies = length(use_energy_inds)
use_energies_ev = energies_kev[use_energy_inds] * 1000
function get_pulses()
    pulses = []
    description = []
    for rn in use_rns
        for i_energy in use_energy_inds
            push!(pulses, pulses_by_rn[rn][1:decimation_factor:end, i_energy])
            push!(description, (;rn=rn,i_energy=i_energy))
        end
    end
    return vcat(pulses...), description
end
ydata, description = get_pulses()
xdata = range(0, step=dt_s, length=length(ydata))
dt = dt_s * 1u"s"

figure()
plot(xdata, ydata * 1e6, ".")
xlabel("times (s)")
ylabel("pulse (uA)")
;
first_pulse_arrival_s = 0.0078667

# define our starting guesses
Rshunt = 360e-3u"mΩ";
Rbias = 1030u"Ω"; Vbias = 0.95u"V"; expected_Vt = Vbias * Rshunt / Rbias |> u"μV";
Tc = 109u"mK"; Tbath = 65.0u"mK"
R0 = 1.0u"mΩ"; Rpara = 1e-3u"mΩ"; Rl = Rshunt + Rpara; Rn = 12.2u"mΩ"; 
n = 3.4; G = 33.3u"pW/K"; C = 0.3u"pJ/K"; L = 50.0u"nH"
target_alpha = 400
target_beta = 1
rit = ShankRIT_from_αβ(target_alpha, target_beta, n, Tc, Tbath, G, R0, Rn)
tes_param = TESParams(n, Tbath, G, C, L, Rl, Rpara, rit)
bt0 = BiasedTES_from_R0(tes_param, R0)
Ic0 = 1u"mA"
rit2 = ModelTES.TwoFluidRIT_from_α(target_alpha, R0, bt0.T0, bt0.I0, Ic0, Rn, Tc)
tes_param2 = TESParams(n, Tbath, G, C, L, Rl, Rpara, rit2)
bt2 = BiasedTES_from_R0(tes_param2, R0)

pulse_arrivals_s = [first_pulse_arrival_s + i * dt_s * n_samples_pulse for i in 0:n_energies - 1]
out = ModelTES.pulses(n_samples_pulse * n_energies, dt_s * 1u"s", bt2, use_energies_ev .* 1u"eV", 
pulse_arrivals_s .* 1u"s", dtsolver=1u"μs")


figure()
plot(ustrip.(u"s", ModelTES.times(out)), ustrip.(u"μA", out.I), label="model guess")
plot(xdata[1:end ÷ n_rns], ydata[1:end ÷ n_rns] * 1e6, label="data")
xlabel("time (s)")
ylabel("pulse (uA)")
legend()
;

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
function make_2fluid_biased_tes(n, Tbath, G, C, L, Rl, Rp, ci, Rn, Ic0, Tc, I0, T0, Vt)
    rit = ModelTES.TwoFluidRIT(ci, Rn * 1u"mΩ", Ic0 * 1u"mA", Tc * 1u"mK")
    p = ModelTES.TESParams(n, Tbath * 1u"mK", G * 1u"pW/K",
      C * 1u"pJ/K", L * 1u"nH", Rl * 1u"mΩ", Rp * 1u"mΩ",  rit)
    ModelTES.BiasedTES(p, I0 * 1u"mA", T0 * 1u"mK", Vt * 1u"mV")
end  
try # allows redefinition of fitfunc more easily
    for m in methods(fitfunc) 
        Base.delete_method(m) 
    end
catch
end
function fitfunc(x, n, Tbath, G, C, L, Rl, Rp, ci, Rn, Ic0, Tc, I0_1, T0_1,
    I0_2, T0_2, I0_3, T0_3, Vt_1, Vt_2, Vt_3, _current_per_arbs)
    bt_1 = make_2fluid_biased_tes(n, Tbath, G, C, L, Rl, Rp, ci, Rn, Ic0, Tc, I0_1, T0_1, Vt_1)
    bt_2 = make_2fluid_biased_tes(n, Tbath, G, C, L, Rl, Rp, ci, Rn, Ic0, Tc, I0_2, T0_2, Vt_2)
    bt_3 = make_2fluid_biased_tes(n, Tbath, G, C, L, Rl, Rp, ci, Rn, Ic0, Tc, I0_3, T0_3, Vt_3)
    y = zeros(n_rns * n_energies * n_samples_pulse)
    for (i, bt) in enumerate([bt_1, bt_2, bt_3][1:n_rns])
        for (j, energy_ev) in enumerate(use_energies_ev)
            out = ModelTES.pulses(n_samples_pulse, dt, bt, energy_ev * 1u"eV",
            [first_pulse_arrival_s * 1u"s"], dtsolver=.1u"μs")
            o = ((i - 1) * n_energies + (j - 1)) * n_samples_pulse
            y[(1:n_samples_pulse) .+ o] = ustrip.(u"A", out.I)
        end
    end
    return y * _current_per_arbs
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
params.C(val=2)
params.ci(val=0.7, vary=true)
params.Ic0(val=1)
params.I0_1(val=pulses_by_rn[use_rns[1]][1,1] * 1e3, vary=false)
params.I0_2(val=n_rns > 1 ? pulses_by_rn[use_rns[2]][1,1] * 1e3 : 0, vary=false)
params.I0_3(val=n_rns > 2 ? pulses_by_rn[use_rns[3]][1,1] * 1e3 : 0, vary=false)
params.T0_1(val=ustrip(u"mK", bt2.T0))
params.T0_2(val=ustrip(u"mK", bt2.T0))
params.T0_3(val=ustrip(u"mK", bt2.T0))
params.Vt_1(val=3.8e-5)
params.Vt_2(val=3.8e-5 * 1.1)
params.Vt_3(val=3.8e-5 * 1.5)
params._current_per_arbs(val=1, vary=false)
out = model(x=xdata, params=params)
;

result = fit(model, params, x=xdata, y=ydata);

figure(figsize=(14, 6))
plot(xdata, ydata * 1e6, label="data")
plot(xdata, result(x=xdata, params=result.params) * 1e6, label="fit")
# plot(xdata, result(x=xdata, params=result.init_params)*1e6, label="guess")
xlabel("time (s)")
ylabel("current (uA)")
legend()
plt.tight_layout()


function ptable(result)
    pretty_table((name = [p.name for p in result.params],
unit = ["", "mK", "pW/K", "pJ/K", "nH", "mΩ", "mΩ", "", "mΩ", "mA", "mK", 
"mA", "mA", "mA", "mK", "mK", "mK", "V", "V", "V", "μA"],
vary = [p.vary for p in result.params],
guess = values(result.init_params),
value = values(result.params),
fit_uncertainty = [p.unc == nothing ? NaN : p.unc for p in result.params],), 
formatters=ft_printf("%5.3g", [3,4,5,6]))
end
ptable(result)
    