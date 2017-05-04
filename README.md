# ModelTES

[![Build Status](https://travis-ci.org/ggggggggg/ModelTES.jl.svg?branch=master)](https://travis-ci.org/ggggggggg/ModelTES.jl)

To install `Pkg.clone("https://github.com/ggggggggg/ModelTES.jl")`

This package is intended to simulate TES microcalorimeter pulses in both the linear and non-linear
regime. It's not really intended for wide use, but feel free to play with it if you want. There may
be no support, but can try opening and issue and hope.

The physics model is based on the article "Transition-Edge Sensors" by
K.D. Irwin and G.C. Hilton in _Cryogenic Particle Detection_ pages 63-150,
Volume 99 of the series _Topics in Applied Physics_ (2005). Find it at
[doi:10.1007/10933596_3](http://doi.org/10.1007/10933596_3).

Here are some examples of what you can do:

```julia
using ModelTES, PyPlot
# Create a high-E TES design
stdTES = ModelTES.highEpix()
# Create a Biased TES from the 48 nanohenry Holmes parameters with 0.2*Rn resistance
tes =  ModelTES.pholmes(48e-9, 0.20)
# Create a Biased TES from the 48 nanohenry Holmes parameters with 0.4*Rn resistance
tes2 = ModelTES.pholmes(48e-9, 0.40)
# Integrate a pulse with 12000 samples, 1e-7 second spacing, 1000 eV energy, 2000 presamples
out = pulse(12000,1e-7, tes, 1000, 2000);
# Integrate a pulse with 12000 samples, 1e-7 second spacing, 1000 eV energy, 2000 presamples from the higher biased version of the same tes
out2 = pulse(12000,1e-7, tes2, 1000, 2000);
# Get all the linear parameters for the irwin hilton model
lintes = IrwinHiltonTES(tes)
# Calculate the noise PSD and its 4 components in the Irwin-Hilton model
f = logspace(0,6,100);
n,n1,n2,n3,n4 = ModelTES.noisePSD(lintes, f);
sampleTime = 1e-6;
nmodel = NoiseModel(tes, sampleTime);
psd = noise_psd(nmodel, f);
covar = model_covariance(nmodel, 100)


# Calculate a stochastic noise 1000 eV pulse with 12000 samples and 2000 presmples
# !! Warning: `stochastic` is not validated.
outstochastic = stochastic(12000,1e-7, tes, 1000,2000);

# calculate the complex impedance of the TES
z = ModelTES.Z(lintes,f)

figure()
title("stochastic pulse and noiseless pulse")
plot(1e3*times(out), 1e6*out.I,".k")
plot(1e3*times(outstochastic), 1e6*outstochastic.I,".r")
xlabel("Time (ms)"); ylabel("Current (\$\\mu\$A)");
figure()
title("pulses with different bias points")
plot(1e3*times(out), 1e6*out.I)
plot(1e3*times(out2), 1e6*out2.I)
xlabel("Time (ms)"); ylabel("Current (\$\\mu\$A)");
figure()
loglog(f,n1,label="TES electrical")
loglog(f,n2,label="Load electrical")
loglog(f,n3,label="Thermal fluctuation")
loglog(f,n4,label="Amplifier")
loglog(f,n,label="Total noise PSD", lw=2)
xlabel("Frequency (Hz)"); ylabel("Noise Power (A\$^2\$/Hz)")
legend(loc="best")
figure()
plot(real.(z), imag.(z))
xlabel("real(z)")
ylabel("imag(z)")
```
