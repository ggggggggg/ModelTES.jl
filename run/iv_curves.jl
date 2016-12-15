using ModelTES#, Plots
using PyPlot
#pyplot() #make Plots use pyplot backend
using PyCall
@pyimport numpy as np

function highEpixHighBeta()
    Tc=0.107; Tbath = 0.065
    R0 = 1.55e-3; Rl = 0.35e-3; Rn = 10.3e-3; Rpara=0.0
    n=3.3; k = 3.50e-9; C = 0.72e-12
    alpha = 92.0; beta = 5; L=300e-9
    model = ModelTES.ShankRIT(alpha, beta, n, Tc, Tbath, k, R0, Rn);
    tes_param = TESParams(n,Tc,Tbath,k,C,L,Rl,Rpara,Rn,model)
    bt = BiasedTES(tes_param, R0)
end



pixel = highEpixHighBeta()

iv_pt = ModelTES.iv_point(pixel.p, pixel.V);
IV_Vs = pixel.V*collect(0.05:0.05:5);
I,T,R,V = ModelTES.iv_curve(pixel.p, IV_Vs);

"calculates T given R and I and the thermal paramters"
function getTgivenRandI(pixel,R,I)
    k = pixel.p.k
    n = pixel.p.n
    Tbath = pixel.p.Tbath
    # I^2*R=k(T^n-Tb^n)
    Tn = Tbath^n-I^2*R/k
    T=Tn^(1/n)
end

using Roots
function getRITgivenV(pixel,V)

    Rl = pixel.p.Rl
    # V=I*(R+Rl)
    fzero((t)->getR0error(t, targetR,p),p.Tc-10*transitionwidth(p), p.Tc+10*transitionwidth(p))


Tgot = getTgivenRandI(pixel,R[1],I[1])
@show Tgot-T[1]



function iv_point_from_T(pixel,T, Iguess=1e-6)
    k = pixel.p.k
    n = pixel.p.n
    Tbath = pixel.p.Tbath
    Rl = pixel.p.Rl
    function thermal(I)
        R = ModelTES.R(I,T,pixel.p)
        I^2*R-k*(T^n-Tbath^n)
    end
    I = fzero(thermal, Iguess, [0, 1])
    R = ModelTES.R(I,T,pixel.p)
    V=I*(R+Rl)
    I,T,R,V
end

iv_point_from_T(pixel,0.106)


# I: current in TES
# T: temperature of TES
# V: thevinen equivalen voltage
# Ibias =  V/Rsh
# Vtes = V-I*Rl
Vtes = V-I*pixel.p.Rl;
# plot(Vtes,I)
# xlabel!("V_{tes} (Volts)")
# ylabel!("Ites (Amps)")
# title!("Tb = $(pixel.p.Tbath) K")

pulserecord = ModelTES.rk8(4000,1e-5,pixel,6000)

Tbaths = 0.05:0.005:0.11
params = deepcopy(pixel.p)
Is = zeros(length(IV_Vs), length(Tbaths));
Ts = zeros(length(IV_Vs), length(Tbaths));
Vs = zeros(length(IV_Vs), length(Tbaths));
Rs = zeros(length(IV_Vs), length(Tbaths));
for (i,Tb) in enumerate(Tbaths)
    params.Tbath=Tb
    I,T,R,V = ModelTES.iv_curve(params, IV_Vs)
    Is[:,i]=I
    Ts[:,i]=T
    Rs[:,i]=R
    Vs[:,i]=V
end
Vtess = Vs-Is*params.Rl;
Ps = Vtess.*Is;
Rsrecalc = Vtess./Is;
Rpercs = Rs/params.Rn;

figure()
plot(Vtess, Is)
xlabel("V_{tes} (Volts)")
ylabel("I_{tes} (Amps)")
legend(loc="best", labels=map(repr, Tbaths))


figure()
plot(Rs,Ts)
xlabel("R_{tes} (Ohms)")
ylabel("T_{tes} (Kelvin)")
ylim(params.Tc-params.RIT.Tw*3, params.Tc+params.RIT.Tw*3)

figure()
plot(Ps,Rs)
xlabel("P_{tes} (Watts)")
ylabel("R_{tes} (Ohms)")
xlim(0,3e-12)
legend(loc="best", labels=map(repr, Tbaths))

# legend(loc="righ", labels=map(repr, Tbaths))
# p = plot(Vtess, Is, leg=false, lab=map(repr, Tbaths'))
# xlabel!("V_{tes} (Volts)")
# ylabel!("I_{tes} (Amps)")
#
# p3 = plot(Vs, Rpercs,  leg=true, lab=map(repr, Tbaths'))
# xlabel!("V_{thevinen} (Volts)")
# ylabel!("R (% Rn)")

# assume R is a function of T only, so R uniquely identifies temperature
# for a given value of R, say 0.4*Rn, we fit P as a function of Tb
"powercurves(Ts, Ps, Ts_chosen)
Ts a 2d array of temperatures from IV curves, where the second index for Tbath, first index is for Vbias
Ps a 2d array of dissipated powers from IV curves, same indicies
Ts_chosen a 1d Vector of temperatures chosen to interpolate power curves on
this may also be used when resistance instead of temperature, relying on the assumption that resistance is only a function of temperature
we know that isn't right so there will be errors assocaited with that"
function powercurves(Ts, Ps, Ts_chosen)
    @show size(Ts)
    @show size(Ps)
    P_interp = zeros(size(Ts,2), length(Ts_chosen))
    @show size(P_interp)
    for (i, T_chosen) in enumerate(Ts_chosen)
        @show i
        PvsTb = [np.interp(T_chosen, Ts[:,i], Ps[:,i]) for i in 1:size(Ps,2)]
        P_interp[:,i]=PvsTb
    end
    P_interp
end


Rpercs_chosen = 0.15:0.25:0.95;
P_interp = powercurves(Rpercs, Ps, Rpercs_chosen)



using LsqFit

Ptofit(Tbaths, T, k, n) =  k*(T^n-Tbaths.^n)
model(x, p) = Ptofit(x,p[1],p[2],p[3])

use_inds = Tbaths.<params.Tc-0.003
Rperc_ind=1
fit = curve_fit(model, Tbaths[use_inds], P_interp[use_inds,Rperc_ind]*1e12, [params.Tc, 3.5e3, 3.3])
Tfit, kfit, nfit = fit.param
Rperc_at_Tfit = Rpercs_chosen[Rperc_ind]




Ts_chosen = [params.Tc-params.RIT.Tw/2, params.Tc, params.Tc+params.RIT.Tw/2]
P_interpT = powercurves(Ts, Ps, Ts_chosen)
T_chosen_ind=1
fitT = curve_fit(model, Tbaths[use_inds], P_interpT[use_inds,T_chosen_ind]*1e12, [params.Tc, 1000, 3])
TfitT, kfitT, nfitT = fitT.param

# p5 = plot(Tbaths, P_interpT*1e12,lab=map(repr, Ts_chosen'),w=2)
# plot!(Tbaths[use_inds], model(Tbaths[use_inds], fit.param), lab="fit",style=:dash,w=2)
# xlabel!("T_{bath} (K)")
# ylabel!("Power required to to get to chosen %Rn (pW)")

j=7 # Tbaths index
TfromTb(Tbath, P, k, n) = (P/k+Tbath^n)^(1/n)
Ps+=1e-30
TrecalcT = [TfromTb(Tbaths[q], Ps[i,q]*1e12, kfitT, nfitT) for i in 1:size(Ps,1), q in 1:size(Ps,2)]
TrecalcR = [TfromTb(Tbaths[q], Ps[i,q]*1e12, kfit, nfit) for i in 1:size(Ps,1), q in 1:size(Ps,2)]



Ts_chosen_2 = [params.Tc-params.RIT.Tw, params.Tc, params.Tc+params.RIT.Tw]
P_interp_correct = powercurves(Ts, Ps, Ts_chosen_2)
figure()
plot(Tbaths, P_interp*1e12, lw=2)
plot(Tbaths[use_inds], model(Tbaths[use_inds], fit.param),"r--",lw=2)
plot(Tbaths, P_interp_correct*1e12,"-.",lw=3, label="correct curve")
xlabel("T_{bath} (Kelvin)")
ylabel("Power to reach resistance (pW)")
legend(loc="best", labels=push!(map(repr, Rpercs_chosen), "fit"))

figure()
plot(Tbaths, P_interpTrecalcR*1e12, label="P_interpTrecalcR")
plot(Tbaths, P_interp[:,Rperc_ind]*1e12, label="PinerpR")
legend(loc="best")
xlabel("Tbath")
ylabel("Power to reach %Rn or Temp (pW)")

P_interpTrecalcT = powercurves(TrecalcT, Ps, [TfitT])
figure()
plot(Tbaths, P_interpTrecalcT*1e12, label="P_interpTrecalcR")
plot(Tbaths, P_interpT[:,T_chosen_ind]*1e12, label="PinerpR")
legend(loc="best")
xlabel("Tbath")
ylabel("Power to reach %Rn or Temp (pW)")

# p6=plot(Ps[:,j]*1e12, Ts[:,j]*1e3,lab="Ts")
# plot!(Ps[:,j]*1e12, TrecalcT[:,j]*1e3, lab="TrecalcT")
# plot!(Ps[:,j]*1e12, TrecalcR[:,j]*1e3, lab="TrecalcR")
# xlabel!("power dissipated on TES (pW)")
# ylabel!("temperature of TES (mK)")
# title!("Tbath = $(Tbaths[j]*1e3) mK")
# print("lastline")
function meshgrid{T}(vx::AbstractVector{T}, vy::AbstractVector{T})
    m, n = length(vy), length(vx)
    vx = reshape(vx, 1, n)
    vy = reshape(vy, m, 1)
    (repmat(vx, m, 1), repmat(vy, 1, n))
end
Is_for_mesh = linspace(0,maximum(Is),90);
Ts_for_mesh = linspace(minimum(Ts), maximum(Ts),100);
Imesh, Tmesh = meshgrid(Is_for_mesh, Ts_for_mesh);
Rdirect = [ModelTES.R(Is_for_mesh[i],Ts_for_mesh[j],params) for i=1:size(Is_for_mesh,1),j=1:size(Ts_for_mesh,1)]

# surface(Imesh,Tmesh,Rdirect ,xlabel="I",ylabel="T",zlabel="R",legend=false)
# surface!(Is, Trecalc, Rs, color=:red)
figure(figsize=(12,8))
surf(Is_for_mesh*1e6,Ts_for_mesh*1e3,1e3*Rdirect')
xlabel("I (uA)")
ylabel("T (mK)")
zlabel("R (mOhm)")
surf(Is*1e6, TrecalcT*1e3, Rs*1e3+0.001, color="cyan")
ax=gca()
ax[:plot](pulserecord.I*1e6, pulserecord.T*1e3, pulserecord.R*1e3, lw=4, color="red")
x1,y1,z1 = pulserecord.I[1]*1e6, pulserecord.T[1]*1e3, pulserecord.R[1]*1e3
x2,y2,z2 = pulserecord.I[end]*1e6, pulserecord.T[end]*1e3, pulserecord.R[end]*1e3
# ax[:quiver](x2,y2,z2,x2-x1,y2-y1,z2-z1)
ax[:plot]([x1,x2],[y1,y2],[z1,z2], color="orange",lw=4)
# ax[:plot]([x1],[y1],[z1], s=[10],color="black")


figure()
Ts_for_thisplot = Ts_chosen
Rs_for_thisplot = [ModelTES.R(Is_for_mesh[i], Ts_for_thisplot[j], params) for i in 1:size(Is_for_mesh,1), j in 1:size(Ts_for_thisplot,1)]
plot(Is_for_mesh, Rs_for_thisplot)
title("R vs T @ fixed I")
xlabel("current (A)")
ylabel("resistance (Ohms)")
legend(labels=map(repr,Ts_for_thisplot), loc="best")
println("lastline")


figure()
Ts_fit_calc = (Ps*1e12/kfit+repmat(Tbaths', size(Ts,1),1).^nfit).^(1/nfit);
surf(Is*1e6, Ts_fit_calc*1e3 ,Rs*1e3)
xlabel("I (uA)")
ylabel("T (mK)")
zlabel("R (mOhm)")

figure()
Ts_fit_calc = (Ps*1e12/kfit+repmat(Tbaths', size(Ts,1),1).^nfit).^(1/nfit);
surf(Is*1e6, Ts*1e3 ,Rs*1e3)
xlabel("I (uA)")
ylabel("T (mK)")
zlabel("R (mOhm)")

figure()
plot(Ts*1e3,  Ts_fit_calc*1e3)
plot()
xlabel("temperature true (mK)")
xlabel("temperature fit_calc (mK)")
legend(loc="best", labels=map(repr, Tbaths),title="Tbath")

figure()
plot(Ts*1e3,  (Ts_fit_calc-Ts)*1e3)
plot([params.Tc-params.RIT.Tw*3, params.Tc+params.RIT.Tw*3]*1e3, [1,1],"k", lw=3)
xlabel("temperature true (mK)")
ylabel("temperature fit_calc-true (mK)")
legend(loc="best", labels=push!(map(repr, Tbaths),"transition"),title="Tbath")
