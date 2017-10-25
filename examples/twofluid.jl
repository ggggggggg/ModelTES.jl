using ModelTES, PyPlot

#following equation 4
# https://doi.org/10.1063/1.4984065 Morgan et al 2017
type TwoFluidRIT <: ModelTES.AbstractRIT
  Ic0::Float64
  cI::Float64
  cR::Float64

end
function ModelTES.R(I,T, RIT::TwoFluidRIT, Tc, Rn)
  a=(1-min(T,Tc)/Tc)^1.5
  b=RIT.cI*RIT.Ic0/I
  c=(1-a*b)
  Rn*RIT.cR*c
end
ModelTES.transitionwidth(RIT::TwoFluidRIT) = 0.001 # I dont think this is used much
tfrit = TwoFluidRIT(0.00742, 0.84, 0.45)
Tbath = 0.055; L = 45e-9; Rl = 0.0; C = 0.25e-12; G=42e-12; n=2.5;
Rn = 0.010; Tc = 0.0714; k=G/n/Tc^(n-1)
R0 = 0.15*Rn; Rpara = 0.0
tes_param = TESParams(n,Tc,Tbath,k,C,L,Rl,Rpara,Rn,tfrit)
bt = BiasedTES(tes_param, R0)

out = pulse(12000,1e-7, bt, 1000, 2000);

plot(times(out), out.I[1]-out.I)
