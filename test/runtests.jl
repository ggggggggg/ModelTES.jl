using ModelTES, Unitful, Base.Test

ModelTES.BiasedTES{ModelTES.ShankRIT}(
ModelTES.TESParams{ModelTES.ShankRIT}(3.15, 0.094, 0.075, 1.48e-9, 1.0e-13, 5.0e-8, 0.00035, 0.0, 0.0082,
ModelTES.ShankRIT(0.0006043070792387806, 0.19880461679172617)),
1.506454430784236e-5, 0.0916286322033873, 2.9978516417800104e-8)

Tc=94.0u"mK"; Tbath = 75.0u"mK"
R0 = 1.64u"mΩ"; Rl = 0.35u"mΩ"; Rn = 8.2u"mΩ"; Rpara=0.0u"mΩ"
n=3.15; G = 27.34u"pW/K"; C = 0.1u"pJ/K"
alpha = 175.0; beta = 2.28; L=50.0u"nH"
model = ShankRIT(alpha, beta, n, Tc, Tbath, G, R0, Rn);
tes_param = TESParams(n,Tbath,G,C,L,Rl,Rpara,model)
I0=1.506454430784236e-2u"mA"
T0= 91.6286322033873u"mK"
V0=2.9978516417800104e-5u"mV"
bt = BiasedTES(tes_param, I0,T0,V0)

Teps = T0*1e-12
@test ModelTES.thermalpower(G, n, T0, T0+Teps) ≈ uconvert(u"W",G*Teps) rtol=1e-2

@show r = ModelTES.R(I0,T0, model)
@show dTout = uconvert(u"K/s",ModelTES.dT(I0, T0, bt.p.G, bt.p.n, bt.p.Tbath, bt.p.C, r))
@show dIout = uconvert(u"A/s",ModelTES.dI(I0,T0, bt.V, bt.p.Rl, bt.p.L, r))
# bt(0,u,du)

out = pulses(1000, 1e-6, bt, [1000u"eV"], [100e-6]; dtsolver=1e-9, method=DifferentialEquations.Tsit5(), abstol=1e-9, reltol=1e-9)
