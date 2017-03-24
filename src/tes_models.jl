using ModelTES

"""Return a BiasedTES object containing parameters for a HOLMES notional TES"""
function pholmes(L=50e-9, RbiasFraction=0.2)
    n = 3.25; Tc = 0.1; Tbath = 0.07
    k = 2.33e-8; C = 0.5e-12
    Rl = 0.3e-3; Rn = 10e-3; Rpara=0.0; R0=Rn*RbiasFraction
    alpha = 200.0; beta = 2.0
    model = ModelTES.ShankRIT(alpha, beta, n, Tc, Tbath, k, R0, Rn);
    tes_param = TESParams(n,Tc,Tbath,k,C,L,Rl,Rpara,Rn,model)
    bt = BiasedTES(tes_param, R0)
end


"""Return a BiasedTES object containing parameters for Kelsey's low-E TES model"""
function lowEpix()
    Tc=0.094; Tbath = 0.075
    R0 = 1.64e-3; Rl = 0.35e-3; Rn = 8.2e-3; Rpara=0.0
    n=3.15; k = 1.48e-9; C = 0.1e-12
    alpha = 175.0; beta = 2.28; L=50e-9
    model = ModelTES.ShankRIT(alpha, beta, n, Tc, Tbath, k, R0, Rn);
    tes_param = TESParams(n,Tc,Tbath,k,C,L,Rl,Rpara,Rn,model)
    bt = BiasedTES(tes_param, R0)
end

"""Return a BiasedTES object containing parameters for Kelsey's high-E TES model"""
function highEpix()
    Tc=0.107; Tbath = 0.065
    R0 = 1.55e-3; Rl = 0.35e-3; Rn = 10.3e-3; Rpara=0.0
    n=3.3; k = 3.50e-9; C = 0.72e-12
    alpha = 92.0; beta = 1.4; L=300e-9
    model = ModelTES.ShankRIT(alpha, beta, n, Tc, Tbath, k, R0, Rn);
    tes_param = TESParams(n,Tc,Tbath,k,C,L,Rl,Rpara,Rn,model)
    bt = BiasedTES(tes_param, R0)
end

"""Return a BiasedTES object containin parameters for the Oct 2016 LCLS-II model"""
function LCLSII(L=52e-9)
    Tc=0.075; Tbath = 0.055
    Rl = 0.34e-3; Rn = 13.5e-3; Rpara=0.0; R0=0.15*Rn
    n=3.5; G=130e-12; k = G/(n*Tc^(n-1)); C = 0.11e-12
    alpha = 120.0; beta = 1.69
    model = ModelTES.ShankRIT(alpha, beta, n, Tc, Tbath, k, R0, Rn);
    tes_param = TESParams(n,Tc,Tbath,k,C,L,Rl,Rpara,Rn,model)
    bt = BiasedTES(tes_param, R0)
end
