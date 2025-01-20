function Weight(proj::Projection,Modeling::AiyagariModelEGM,Solution::AiyagariSolution)
############################################################

@unpack params,ns,states = Modeling
@unpack β,α,τ,  = params

@unpack R,w,TT,K,Ltot = Solution

@unpack xsis, Nbin,ytype,Sp = proj
@unpack xsiu0,xsiu1,xsiuE,xsiv0,xsiv1,xsiy,xsil = xsis
abp = proj.abp./Sp #per capita
lp = proj.lp./Sp #per capita
PI = proj.Matab #transition matrix
uc = proj.uc #per capita
ucc = proj.ucc #per capita

fvp = Modeling.Functions.fvp
fvpp = Modeling.Functions.fvpp

FL = (1-α)*K^α*Ltot^-α

xsiu1t = xsiu1./lp
xsiv1t = xsiv1 ./((1-τ)*w*xsiy.*(states[ytype].^(1-τ)).*(lp.^(-τ)) )
xsiv0t = xsiv0 ./((1-τ)*w*xsiy.*(states[ytype].^(1-τ)).*(lp.^(-τ)) )


#######################
# Constructing Matrices
#######################

Sp[findall(x->x<=10^-8,Sp)] .=10^-8
Diagp = diagm(Sp)
iDiagp = diagm(1.0 ./Sp)

PIt = transpose(PI)
DS = diagm(Sp)

PIb = Diagp*transpose(PI)*iDiagp #per capita
#PIb= iDiagp*PI*Diagp
Id = sparse(I,Nbin,Nbin)

P = sparse(I,proj.Nbin,proj.Nbin);
for i=proj.indcc
    P[i,i]=0
end

Pc = Id - P
onev = ones(Nbin,1)

M1 = inv(diagm(xsiv1t.*fvpp(lp) + τ*xsiu1t.*uc))
M0 = -M1*diagm(xsiv0t.*fvp(lp))
V0 = (FL*M1*(Sp.*(states[ytype])))./((1-τ)*w*xsiy.*(states[ytype].^(-τ+1)).*lp.^(-τ))

M0h = diagm(xsiu0.*uc)
M1h = -diagm(xsiuE.*ucc)*(Id-R*PI)
M2h = (1-τ)*w*diagm(xsiy.*((states[ytype].*lp).^(1-τ)).*xsiu1t.*ucc)

M2 = Id - M2h*M1
M3 = inv(M2)*(M0h + M2h*M0)
M4 = inv(M2)*M1h
V1 = inv(M2)*(M2h*V0 - Sp)

R5t = - inv((Id - P) +P*(Id - β*R*PIb)*M4)*P*(Id - β*R*PIb)
M5 = R5t*M3
V2 = R5t*V1

C1 =  transpose(abp)*(M4*V2+V1) + transpose(xsiuE.*uc)*PI*V2
C1 = C1[1]
L1 = (transpose(abp)*(M3+M4*M5) + transpose(xsiuE.*uc)*PI*M5)/C1

M6 = M3+M4*(M5 - V2*L1) - V1*L1
M6h = M0 +M1*M6-V0*L1


#######################
# Constraints
#######################

L1h = UGp.+ L1
L2 = transpose(log.(states[ytype].*lp).*xsiy.*(states[ytype].^(1-τ).*lp.^(1-τ)))*M6 + transpose((onev+(1-τ)*log.(states[ytype].*lp)).*xsiy.*(states[ytype].^(1-τ) .*lp.^(1-τ)).*xsiu1t.*uc)*M6h
L3 = transpose(xsiy.*states[ytype].^(1-τ).*lp.^(1-τ))*M6 + (1-τ)*transpose(xsiy.*states[ytype].^(1-τ).*lp.^(1-τ).*xsiu1t.*uc)*M6h
L4 = transpose(onev)*M6
L5 = transpose(onev)


#######################
# Weights
#######################

M7 = 1.0

iDy = diagm(1.0 ./ (Sp))
#iDy = 1.0

if UG == 0
    if TT == 0.0
        Vtc = [L2*DS*M7;L3*DS*M7;L5*DS*M7]
        M8 = [L2*DS*M7;L3*DS*M7; L5*DS*M7]*[iDy*transpose(L2*DS*M7) iDy*transpose(L3*DS*M7) iDy*transpose(L5*DS*M7) ]
        V8 = [-(L2*DS*M7*ones(Nbin,1))[1],-(L3*DS*M7*ones(Nbin,1))[1], -(L5*DS*M7*ones(Nbin,1) .-1)[1]]
        muvc = M8\V8
        weights = ones(Nbin,1) + muvc[1]*iDy*transpose(L2*DS*M7)+ muvc[2]*iDy*transpose(L3*DS*M7) + muvc[3]*iDy*transpose(L5*DS*M7)
    else
        Vtc = [L2*DS*M7;L3*DS*M7;L4*DS*M7;L5*DS*M7]
        M8 = [L2*DS*M7;L3*DS*M7;L4*DS*M7;L5*DS*M7]*[iDy*transpose(L2*DS*M7) iDy*transpose(L3*DS*M7) iDy*transpose(L4*DS*M7) iDy*transpose(L5*DS*M7)]
        V8 = [-(L2*DS*M7*ones(Nbin,1))[1],-(L3*DS*M7*ones(Nbin,1))[1], -(L4*DS*M7*ones(Nbin,1))[1], -(L5*DS*M7*ones(Nbin,1) .-1)[1]]
        muvc = M8\V8
        weights = ones(Nbin,1) + muvc[1]*iDy*transpose(L2*DS*M7) + muvc[2]*iDy*transpose(L3*DS*M7) + muvc[3]*iDy*transpose(L4*DS*M7) + muvc[4]*iDy*transpose(L5*DS*M7)
    end
else
    if TT == 0.0
        Vtc = [L1h*DS*M7;L2*DS*M7;L3*DS*M7; L5*DS*M7]
        M8 = [L1h*DS*M7;L2*DS*M7;L3*DS*M7; L5*DS*M7]*[iDy*transpose(L1h*DS*M7) iDy*transpose(L2*DS*M7) iDy*transpose(L3*DS*M7) iDy*transpose(L5*DS*M7)]
        V8 = [-(L1h*DS*M7*ones(Nbin,1))[1],-(L2*DS*M7*ones(Nbin,1))[1],-(L3*DS*M7*ones(Nbin,1))[1], -(L5*DS*M7*ones(Nbin,1) .-1)[1]]
        muvc = M8\V8
        weights = ones(Nbin,1) +muvc[1]*iDy*transpose(L1h*DS*M7)+ muvc[2]*iDy*transpose(L2*DS*M7)+ muvc[3]*iDy*transpose(L3*DS*M7) + muvc[4]*iDy*transpose(L5*DS*M7)
    else
        Vtc = [L1h*DS*M7;L2*DS*M7;L3*DS*M7;L4*DS*M7 ; L5*DS*M7]
        M8 = [L1h*DS*M7;L2*DS*M7;L3*DS*M7;L4*DS*M7 ; L5*DS*M7]*[ iDy*transpose(L1h*DS*M7) iDy*transpose(L2*DS*M7) iDy*transpose(L3*DS*M7) iDy*transpose(L4*DS*M7) iDy*transpose(L5*DS*M7) ]
        V8 = [-(L1h*DS*M7*ones(Nbin,1))[1],-(L2*DS*M7*ones(Nbin,1))[1],-(L3*DS*M7*ones(Nbin,1))[1], -(L4*DS*M7*ones(Nbin,1))[1], -(L5*DS*M7*ones(Nbin,1) .-1)[1]]
        muvc = M8\V8
        weights = ones(Nbin,1) + muvc[1]*iDy*transpose(L1h*DS*M7)+ muvc[2]*iDy*transpose(L2*DS*M7) + muvc[3]*iDy*transpose(L3*DS*M7) + muvc[4]*iDy*transpose(L4*DS*M7) + muvc[5]*iDy*transpose(L5*DS*M7)
    end
end


if positive == true
    x = findall(x->x<=0,weights)
    wmin = 0.01
    count=1

    while ((!isempty(x))&(count<=10))
        Sx = length(x)
        Lv = zeros(Sx,Nbin)
        Cv = zeros(Nbin,Sx)
        for i=1:Sx
            Lv[i,x[i]] = 1
        end

        for i=1:Sx
            Cv[:,i] .= iDy*(Lv[i,:])
        end


        if UG == 0
            if TT == 0.0
                M8 = [L2*DS*M7;L3*DS*M7;L5*DS*M7;Lv]*[iDy*transpose(L2*DS*M7) iDy*transpose(L3*DS*M7) iDy*transpose(L5*DS*M7) Cv]
                V8x =  wmin*ones(Sx,1)- Lv*ones(Nbin,1)
                V8 = [[-(L2*DS*M7*ones(Nbin,1))[1],-(L3*DS*M7*ones(Nbin,1))[1],-(L5*DS*M7*ones(Nbin,1) .- 1)[1]];V8x]
                muvc = M8\V8
                weightx = muvc[4:end]
                weights = ones(Nbin,1) + muvc[1]*iDy*transpose(L2*DS*M7)+ muvc[2]*iDy*transpose(L3*DS*M7) + muvc[3]*iDy*transpose(L5*DS*M7) + Cv*weightx
            else
                M8 = [L2*DS*M7;L3*DS*M7;L4*DS*M7;L5*DS*M7;Lv]*[iDy*transpose(L2*DS*M7) iDy*transpose(L3*DS*M7) iDy*transpose(L4*DS*M7) iDy*transpose(L5*DS*M7) Cv]
                V8x =  wmin*ones(Sx,1)- Lv*ones(Nbin,1)
                V8 = [[-(L2*DS*M7*ones(Nbin,1))[1],-(L3*DS*M7*ones(Nbin,1))[1], -(L4*DS*M7*ones(Nbin,1))[1], -(L5*DS*M7*ones(Nbin,1) .- 1)[1]];V8x]
                muvc = M8\V8
                weightx = muvc[5:end]
                weights = ones(Nbin,1) + muvc[1]*iDy*transpose(L2*DS*M7) + muvc[2]*iDy*transpose(L3*DS*M7) + muvc[3]*iDy*transpose(L4*DS*M7) +  muvc[4]*iDy*transpose(L5*DS*M7) + Cv*weightx
            end
        else
            if TT == 0.0
                M8 = [L1h*DS*M7;L2*DS*M7;L3*DS*M7;L5*DS*M7;Lv]*[iDy*transpose(L1h*DS*M7) iDy*transpose(L2*DS*M7) iDy*transpose(L3*DS*M7) iDy*transpose(L5*DS*M7) Cv]
                V8x =  wmin*ones(Sx,1)- Lv*ones(Nbin,1)
                V8 = [[-(L1h*DS*M7*ones(Nbin,1))[1],-(L2*DS*M7*ones(Nbin,1))[1],-(L3*DS*M7*ones(Nbin,1))[1],-(L5*DS*M7*ones(Nbin,1) .- 1)[1]];V8x]
                muvc = M8\V8
                weightx = muvc[5:end]
                weights = ones(Nbin,1) +muvc[1]*iDy*transpose(L1h*DS*M7)+ muvc[2]*iDy*transpose(L2*DS*M7)+ muvc[3]*iDy*transpose(L3*DS*M7) + muvc[4]*iDy*transpose(L5*DS*M7) + Cv*weightx
            else
                M8 = [L1h*DS*M7;L2*DS*M7;L3*DS*M7;L4*DS*M7; L5*DS*M7; Lv]*[ iDy*transpose(L1h*DS*M7) iDy*transpose(L2*DS*M7) iDy*transpose(L3*DS*M7) iDy*transpose(L4*DS*M7) iDy*transpose(L5*DS*M7) Cv]
                V8x =  wmin*ones(Sx,1)- Lv*ones(Nbin,1)
                V8 = [[-(L1h*DS*M7*ones(Nbin,1))[1],-(L2*DS*M7*ones(Nbin,1))[1],-(L3*DS*M7*ones(Nbin,1))[1], -(L4*DS*M7*ones(Nbin,1))[1], -(L5*DS*M7*ones(Nbin,1) .- 1)[1]];V8x]
                muvc = M8\V8
                weightx = muvc[6:end]
                weights = ones(Nbin,1) + muvc[1]*iDy*transpose(L1h*DS*M7)+ muvc[2]*iDy*transpose(L2*DS*M7) + muvc[3]*iDy*transpose(L3*DS*M7) + muvc[4]*iDy*transpose(L4*DS*M7) +  muvc[5]*iDy*transpose(L5*DS*M7) + Cv*weightx
            end
        end

        #x = findall(x->x<=0,weights)
        y = findall(x->x<=0,weights)
        if !isempty(y)
                x = unique((vcat(x,y)))
        end
        count = count+1

    end
else
    weights = weights
end


#######################
# Solution
#######################

weightb = DS*M7*weights
mu =  - L1*weightb
mu = mu[1]
mut = sum(UGp.*weightb)
ψb = M6*weightb
lambdacb =  (M5 - V2*L1)*weightb
lambdalb = M0*weightb + M1*ψb + mu*V0
lambdacbt = PI*lambdacb

weight2 = ones(1,Nbin) #pareto weights with full instrments

CC = zeros(8)

CC[1] = maximum(abs.(-ψb  + weightb.*xsiu0.*uc - (lambdacb .*xsiuE- R.*(PI*lambdacb).*xsiuE - (1-τ)*w.*lambdalb.*xsiy.*xsiu1t.*(states[ytype].^(1-τ).*lp.^(1-τ))).*ucc - mu.*Sp))
CC[2] = maximum(abs.(P*(ψb -  β*R*PIb*ψb)))
CC[3] = maximum(abs.(Pc*lambdacb))
CC[4] = maximum(abs.( - (weightb.*xsiv0t.*fvp(lp) + xsiv1t.*lambdalb.*fvpp(lp))  + ψb + (-τ)*xsiu1t.*uc.*lambdalb  + (mu*FL/((1-τ)*w))*Sp./(xsiy.*states[ytype].^(-τ).*lp.^(-τ)) ))
CC[5] = maximum(abs.(transpose(abp)*ψb + transpose(xsiuE.*uc)*PI*lambdacb))
CC[6] = maximum(abs.(transpose(xsiy.*states[ytype].^(1-τ).*lp.^(1-τ))*ψb + (1-τ)*transpose(xsiy.*uc.*xsiu1t.*states[ytype].^(1-τ).*lp.^(1-τ))*lambdalb))
CC[7] = maximum(abs.(transpose(log.(states[ytype].*lp).*xsiy.*(states[ytype].^(1-τ).*lp.^(1-τ)))*ψb + (transpose( (onev+(1-τ)*log.(states[ytype].*lp)).*xsiy.*(states[ytype].^(1-τ) .*lp.^(1-τ)).*xsiu1t.*uc)*lambdalb)))
CC[8] = sum(ψb)


res = Planner(
vec(lambdacb),
vec(lambdalb),
vec(lambdacbt),
vec(ψb),
vec(weight2),
vec(weights),
vec(weightb),
mu)

J1 = L1h*DS*M7
J2 = L2*DS*M7
J3 = L3*DS*M7
J4 = L4*DS*M7
J5 = L5*DS*M7

#VS = 1/(M7'*(Sp.*uc)./disty) #pareto weight with full instrments

return res,CC, J1, J2, J3, J4, J5
end
