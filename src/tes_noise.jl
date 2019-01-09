using Polynomials
using ARMA

"""`noisePSD(tes, freq, [SI_amp])`

Returns (noise,A,B,C,D) where `noise` is a noise power spectral
density in A^2 per Hz, at each input frequency (Hz) given by `freq`.
A,B,C,D are the components of the noise, such that A+B+C+D=noise.
In order, they are TES Johnson noise; load resistor Johnson noise;
thermal fluctuation noise; and amplifier noise.

You can input the amplifier current noise `SI_amp`, in A^2/Hz, or you
can take the default value."""
function noisePSD(tes::IrwinHiltonTES, freq::AbstractVector, SI_amp=5e-22)
    # This term F goes from 0 to one and depends on whether the thermal conductivity is ballistic
    # or diffusive, hardcoded as 1 for now
    F  = 1
    SP_TFN = 4kb*tes.T0^2*tes.G0*F
    SV_TES = 4kb*tes.T0*tes.R0*(1+2*tes.beta) # TES voltage noise
    SV_L   = 4kb*tes.T0*tes.Rl  # Load voltage noise

    omega = 2pi*freq  # Radians / sec
    sIomeg = (1-tes.tauplus/tes.taucc)*(1-tes.tauminus/tes.taucc)./((1.+im*omega*tes.tauplus).*(1.+im*omega*tes.tauminus)) /(tes.I0*tes.R0*(2+tes.beta))
    sIomeg2 = abs2.(sIomeg)

    Inoise_TES = SV_TES*tes.I0^2/tes.loopgain^2 * (1 .+(tes.tauthermal*omega).^2) .* sIomeg2
    Inoise_load = SV_L*tes.I0^2*(tes.loopgain-1)^2/tes.loopgain^2 * (1 .+(tes.taucc*omega).^2) .* sIomeg2
    Inoise_TFN = SP_TFN*sIomeg2
    Inoise_amp = SI_amp

    Inoise = Inoise_TFN+Inoise_amp+Inoise_TES+Inoise_load
    Inoise, Inoise_TES, Inoise_load, Inoise_TFN, Inoise_amp .+ zeros(Float64, length(Inoise))
end

noisePSD(tes::BiasedTES, freq::AbstractVector, SI_amp::Float64=5e-22) =
    noisePSD(IrwinHiltonTES(tes), freq, SI_amp)


"""`NoiseModel(tes, sampleTime, [SI_amp])`

Returns noise of an IrwinHiltonTES modeled as an ARMA(2,2) process, specifically as an
ARMAModel object from package ARMA.jl"""
function NoiseModel(tes::IrwinHiltonTES, sampleTime::Float64, SI_amp=5e-22)
    # This term F goes from 0 to one and depends on whether the thermal conductivity is ballistic
    # or diffusive, hardcoded as 1 for now
    F  = 1
    SP_TFN = 4kb*tes.T0^2*tes.G0*F
    SV_TES = 4kb*tes.T0*tes.R0*(1+2*tes.beta) # TES voltage noise
    SV_L   = 4kb*tes.T0*tes.Rl  # Load voltage noise

    A = SI_amp
    B = abs2((1-tes.tauplus/tes.taucc)*(1-tes.tauminus/tes.taucc) /(tes.I0*tes.R0*(2+tes.beta)))
    C = SP_TFN
    D = SV_TES*tes.I0^2/tes.loopgain^2
    E = SV_L*tes.I0^2*(tes.loopgain-1)^2/tes.loopgain^2

    # (u,v) are the solutions to the quadratic equation ax^2+bx+c=0
    a = A * (tes.tauplus*tes.tauminus)^2
    b = -(A*(tes.tauplus^2+tes.tauminus^2)+B*(D*tes.tauthermal^2+E*tes.taucc^2))
    c = A+B*(C+D+E)
    discr = sqrt(b^2-4a*c)
    u2 = (-b-discr)/(2a)
    v2 = (-b+discr)/(2a)

    u,v = sqrt(u2), sqrt(v2)
    K = sqrt(a)
    theta = [(sampleTime*u*.5+1)*(sampleTime*v*0.5+1),
            sampleTime^2*u*v*.5-2,
            (sampleTime*u*.5-1)*(sampleTime*v*.5-1)]
    phi = [(2tes.tauplus/sampleTime+1)*(2tes.tauminus/sampleTime+1),
            2-8tes.tauplus*tes.tauminus/sampleTime^2,
            (1-2tes.tauplus/sampleTime)*(1-2tes.tauminus/sampleTime)]

    # Rescale so that phi has unit coefficients for the z^0 term; theta[1] is sigma.
    sigma = theta[1]/phi[1]*4K/(sampleTime^2)
    sigma /= sqrt(2sampleTime)
    theta .*= sigma/theta[1]
    phi ./= phi[1]

    ARMAModel(real(theta), real(phi))
end

NoiseModel(tes::BiasedTES, sampleTime::Float64, SI_amp::Float64=5e-22) =
    NoiseModel(IrwinHiltonTES(tes), sampleTime, SI_amp)
