using ModelTES
using Base.Test

Tc=0.107; Tbath = 0.065
R0 = 1.55e-3; Rl = 0.35e-3; Rn = 10.3e-3; Rpara=0.0
n=3.3; k = 3.50e-9; C = 0.72e-12
alpha = 92.0;β = 1.4; L=300e-9
model = ModelTES.ShankRIT(alpha, β, n, Tc, Tbath, k, R0, Rn);
tes_param = TESParams(n,Tc,Tbath,k,C,L,Rl,Rpara,Rn,model)
bt = BiasedTES(tes_param, R0)
ihtes = IrwinHiltonTES(bt)

@test Tbath == ihtes.Tbath
@test isapprox(R0, ihtes.R0; rtol=1e-5)
@test Rl == ihtes.Rl
@test C == ihtes.C0
@test isapprox(alpha, ihtes.alpha; rtol=1e-5)
@test isapprox(β, ihtes.beta; rtol=1e-5)
@test L == ihtes.L

f = logspace(0,6,100);
n,n1,n2,n3,n4 = ModelTES.noisePSD(ihtes, f);
z = ModelTES.Z(ihtes,f)
zcircuit = ModelTES.Zcircuit(ihtes,f)

zr, zi = real.(z), imag.(z)

# Check issue #4: is the normalization correct on PSD and autocorrelation?
# Estimate the total power by a weird log-space trapezoid rule.
amplifier_noise=1e-35
frequencies = collect(logspace(0,6,1501));
psdRef = ModelTES.noisePSD(bt, frequencies, amplifier_noise)[1];
totalpower = psdRef[1]*frequencies[1];
for i=2:length(frequencies)
    totalpower += (frequencies[i]-frequencies[i-1])*0.5*(psdRef[i]+psdRef[i-1]);
end

sampleTime = 1e-8
modelRef = NoiseModel(bt, sampleTime, amplifier_noise)
for sampleTime in logspace(-7.5,-6,4)
    model = NoiseModel(bt, sampleTime, amplifier_noise)
    @test isapprox(model.covarIV[1], modelRef.covarIV[1], rtol=0.01)
    @test isapprox(model.covarIV[1], totalpower, rtol=0.2)
end
