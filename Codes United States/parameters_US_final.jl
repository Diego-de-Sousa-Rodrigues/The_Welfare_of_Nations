##############################
#   Defining parameters
##############################

#############Parameters
α = 0.36 #capital-share
δ = 0.025 #depreciation in quarterly terms
γ = 1.8 #concavity of the utility function
ϕ = .5
tk = 0.36 #taxes on capital
tc = 0.05 #taxes on consumption
κ = .85 #initial guess
τ = 0.16 #be careful, it is 1-τ in Heathcote and κ is λ !
Tr_init = 0 # should try to reach 8%

#############Variables in the steady-state
KsY = 2.7*4 #capital to output ratio
β = (1+KsY^(-1)*α-δ)^-1 #discount factor
Rt = 1/β #before-tax gross interest rate
wt = (1-α)*( (Rt - (1-δ))/α )^(α/(α-1))
R = (1 + (1 - tk)*(Rt-1))*(1+tc)./(1+tc)
w = (κ*wt^(1-τ))./(1+tc)
χ = ((1/3.5)^(1/ϕ))/w

TT = Tr_init/(1+tc)


#############Idiosyncratic process - United States. Those values were obtained in the paper "Understanding cross country differences in health status and expenditures"
pers = 0.96184
sigma2=sqrt(0.03845)
vary = (sigma2^2/(1 - pers^2))

rho    = pers^0.25; #persistence component of income
sige = ((sigma2^2/(1 + rho^2 + rho^4 + rho^6 ))^(0.5))


#########################################################
 # Functions
#########################################################

#############Utility function and derivatives
(γ ==1.0) ? fu(c)=log.(c) : fu(c)=(c.^(1-γ))./(1-γ)
(γ ==1.0) ? fuinv(X) = exp.(X) :  fuinv(X) = ((1-γ).* ( X ) ).^(1/(1-γ))
fup(c) = c.^-γ
fupp(c) = -γ*c.^(-γ-1)
fv(l) = (1/χ)*(l.^(1+1/ϕ))./(1+1/ϕ) # function v
fvp(l) = (1/χ)*l.^(1/ϕ)
fvpp(l) = (1/χ) *(1/ϕ)*l.^(1/ϕ-1)

Functions = Aiyagarifunction(fu,fuinv,fup,fupp,fv,fvp,fvpp)
ns     = 10
na     = 100
amin = 1e-9
amax = 1000.0
curv = 4
aGrid  = (grid_fun(amin,amax,na,curv))./(1+tc)


#############Idiosyncratic process - Tauchen
mc     = tauchen(ns,rho,sige)
endow  = exp.(mc.state_values)
Trans = collect(mc.p')
Trans[findall(x->x<=5*10^-5,Trans)].=0
for i =1:ns
    Trans[i,i] = Trans[i,i]+ (1 - sum(Trans,dims=1)[i])
end
disty  = (Trans^100000)[:,1] #stationary distribution for the idosyncratic process
Transv = vec(Trans)


#############Renormalization
states = endow./dot(disty,endow)
yvect = vec(kron(states,ones(na,1)));


#############Results
params = AiyagariParametersEGM(β,α,δ,γ,χ,ϕ,κ,τ,tc,tk,TT)
pol0   = repeat(aGrid,ns)
Modeling  = AiyagariModelEGM(params, Functions, aGrid,na,states,ns,Transv,disty,yvect)
