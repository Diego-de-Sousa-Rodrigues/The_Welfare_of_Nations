#########################################################
 # Projection for the Truncation Method
#########################################################


function  convertBasisE10(vE::Array{Int64,1}, E::Int64, N::Int64)
    #converts the vector vE into an integer p
    #vE is of length N, containing integers between 1 and E
    #p is the decimal represention of v (which is in basis E)
    #convention: vE(1) represents the idiosyncratic state in the present period and vE(N) represents the past idosyncratic state N-1 periods before

    if (length(vE)!=N)
        println("incompatible length for vE");
    end

    if minimum(vE)<1
        println("all elements of vE should be >= 1");
    end

    if maximum(vE)>E
        println("all elements of vE should be <= E");
    end

    p = 0;
    for k=1:N
        p = p+(vE[k]-1)*E^(N-k);
    end

    p=p+1; #otherwise between 0 and (E+1)^N - 1
    return p
end


function  convertBasis10E(p::Int64, E::Int64, N::Int64)
    #converts the integer p into a vector vE of length N, containing integers between 1 and E+1 representing p in basis E
    #convention: vE(1) represents the idiosyncratic  state in the present period and vE(N) represents the past idosyncratic state N-1 periods before

    vE=Array{Int64,1}

    if ((p<1)|(p>E^N))
        println("incompatible integer p");
    end

    ptemp = p-1;
    vE = zeros(Int64,1,N);r=0;
    for k=1:N
        r = mod(ptemp, E);
        vE[N-k+1] = round(Int,1+r);
        ptemp = floor(Int,(ptemp-r)/(E));
    end
    return vE
end


#########################################################
 # Main Function
#########################################################


function Project_plan(
    N::Integer, #length of the truncation
    Model::AiyagariModelEGM,
    Solution::AiyagariSolution)

        @unpack params,aGrid,na,ns,Transv,states,disty = Modeling
        @unpack γ,β,χ,ϕ,τ = params
        @unpack ga,gl,gc,dist,Mat,resE,R,w,TT = Solution

        fu = Modeling.Functions.fu #utility of consumption
        fup = Modeling.Functions.fup #marginal utility of consumption
        fv = Modeling.Functions.fv #disutility of labor
        fvp = Modeling.Functions.fvp #marginal disutility of labor

        Trans   = reshape(Transv,ns,ns)
        ytype = kron(collect(StepRange(1, Int8(1), ns)),ones(Int64,na)) #on the initial grid define the types of agents in this economy

        upc  = fup(gc)
        upcc = fupp(gc)

        lhsvi = (states[ytype].*gl).^(1-τ) # term entering the individual budget constraint
        v0 = fv.(gl)
        v1 = fvp.(gl)

        Nidio = ns #this is to update the sizes of the bin according to the function we had defined previously
        Size = zeros(Float64,Nidio^N) #size of the bin
        compt=0
        Ntoi = zeros(Int,Nidio^N)
        for i=1:Nidio^N
            vE =  convertBasis10E(i, Nidio, N)
            s0 = disty[vE[N]] #initial number by productity types
            Siz = s0
            for j=1:N-1
                Siz = Siz*Trans[vE[N-j],vE[N+1-j]]
            end
            if Siz>0
                Size[compt+1]=Siz
                Ntoi[compt+1]=i
                compt+=1
            end
        end
        @show Nbin = compt

        ytypeb = zeros(Int64,Nbin)

        for k=1:Nbin
            i = Ntoi[k] #number of the history
            vE =  convertBasis10E(i, Nidio, N) #this is to recover the vector of histories in which the Size is greater than zero
            ytypeb[k] = mod(vE[1] - 1 ,ns)+1 #productivity type between 1 and ns
        end

        CCv = zeros(Float64,na) #to compute the fraction of credit-constrained bin
        CCv[1] = 1
        CCv = repeat(CCv,Nidio)  #vector of 1 only for the lowest wealth bin for each Nidio

        dGridl_t = repeat(aGrid,ns)
        conso   = zeros(Float64,Nbin) #consumption in the bin
        at   = zeros(Float64,Nbin) #initial of period assets
        ae   = zeros(Float64,Nbin) #this represents the policy function for total assets
        ae2   = zeros(Float64,Nbin) #this represents the policy function for total assets
        uc   = zeros(Float64,Nbin) #marginal utility of the bin (per capita)
        ucc   = zeros(Float64,Nbin) #second derivative of the utility of the bin (per capita)
        Euc  = zeros(Float64,Nbin) #expected utility of the bin
        ls   = zeros(Float64,Nbin) #labor supply of the bin
        lhsv   = zeros(Float64,Nbin) #term entering the budget constraint for the bin
        ul   = zeros(Float64,Nbin) #utility level of conso of the bin
        xsiU = zeros(Float64,Nbin)
        xsiUc= zeros(Float64,Nbin)
        inc  = zeros(Float64,Nbin) #income of the agents
        resp = zeros(Float64,Nbin) #fraction of agents credit constrained
        CC = zeros(Float64,Nbin) #share of credit constrained agents
        X = zeros(Float64,Nbin) #consumption in per capita terms
        v0b = zeros(Float64,Nbin) #labor disutility of the bin (per capita)
        v1b = zeros(Float64,Nbin) #marginal labor disutility of the bin (per capita)
        Transu = spzeros(Nbin,Nbin); #general transition matrix agents
        Sp = zeros(Float64,Nbin) #size of the bin

        Li= zeros(Float64,Nbin)
        Ri= zeros(Float64,Nbin)
        Uc= zeros(Float64,Nbin)
        Vl= zeros(Float64,Nbin)

        for k=1:Nbin
            i = Ntoi[k] #number of the history
            vE =  convertBasis10E(i, Nidio, N) #this is to recover the vector of histories in which the Size is greater than zero
            yy = mod(vE[1] - 1 ,ns)+1 #productivity type between 1 and ns
            id0 = 1+ (vE[N]-1)*na  #index of the initial distribution in dist
            if0 =  (vE[N])*na  #final index of the initial distribution in dist
            DD = 0 .*dist
            DD[id0:if0] .= dist[id0:if0] #this is the initial distribution of agents with the beginning of history 4 for example, considering the case in which we have k = 18 and in this case we have i = 28 and vE= [5 4]
            for j=1:N-1
                Mat_trans = 0 .*Mat
                id0 = 1+ (vE[N+1-j]-1)*na
                if0 =  (vE[N+1-j])*na
                id1 = 1+ (vE[N-j]-1)*na
                if1 =  (vE[N-j])*na
                Mat_trans[id1:if1,id0:if0] = Mat[id1:if1,id0:if0]
                DD = Mat_trans*DD
         end
         Sp[k]         = sum(DD) #this is to check whether the size we obtained above is the same as this one
         CC[k]         = sum(CCv.*(Mat*DD))/Sp[k] #fraction of CC agent in k
         ls[k]         = sum(DD.*gl) #not per capita
         lhsv[k]       = sum(DD.*lhsvi) #not per capita
         #CC[k]        = sum(resE.*DD)/Sp[k] #fraction of CC agent in k
         conso[k]      = sum(DD.*gc) # not percapita
         at[k]         = sum(dGridl_t.*DD) #beginnning of period
         ae2[k]        = sum(dGridl_t.*(Mat*DD)) #total assets in the bin k
         ae[k]         = at[k]*R + w*lhsv[k] + Sp[k]*TT - conso[k] #this represents the total assets

         uc[k]         = sum(upc.*DD)/ Sp[k]
         ucc[k]        = sum(upcc.*DD)/ Sp[k]
         #uc[k]         = sum(upc.*DD)
         #ucc[k]        = sum(upcc.*DD)

         Euc[k]        = sum(upc.*(Mat*DD)) #following utility
         X[k]          = conso[k]/Sp[k] #consumption per capita
         #X[k]          = conso[k] #conso per capita

         #ul[k]         = sum(DD.*fu(gc))
         #v0b[k]        = sum(DD.*v0)
         #v1b[k]        = sum(DD.*v1)
         #resp[k]       = sum(DD.*resE)
         ul[k]         = sum(DD.*fu(gc))/ Sp[k]
         v0b[k]        = sum(DD.*v0)/ Sp[k]
         v1b[k]        = sum(DD.*v1)/ Sp[k]
         resp[k]       = sum(DD.*resE)/ Sp[k]


         Li[k] = (w*lhsv[k])/Sp[k] #labor income after taxes per capita
         Ri[k] = (at[k]*(R-1))/Sp[k] #capital income after taxes per capita
         Uc[k] = sum(DD.*fu(gc))/Sp[k] #utility of consumption
         Vl[k] = sum(DD.*fv(gl))/Sp[k] #disutility of labor
     end


     #############Transition Matrix for agents and wealth
     #constructing the transition matrices on the projection
     Transp       = spzeros(Nbin,Nbin); # general transition matrix AGENTS

     for k  = 1:Nbin
         i  = Ntoi[k] #number of the history
         v  = convertBasis10E(i,Nidio,N)
         statei   = v[1]
         vf = v                     #to initialize
         vf[2:end] = v[1:end-1]     #translated vector
         for statej=1:Nidio
             if Trans[statej,statei]>0 #if possible continuation
                 vf[1] = statej
                 iff = convertBasisE10(vec(vf),Nidio,N)
                 j = findall(x->x==iff,Ntoi) #index of the history
                 if isempty(j)
                 else j=j[1]
                     Transp[j,k] = Trans[statej,statei] #from k to j per capita
                     #Transu[j,k] = Transu[j,k]/uc[j]
                 end
             end
         end
     end


     #############Define the number of credit constrainted bins
     DDD = reshape(dist,na,ns) #disitribution of wealth
     Sh = sum(DDD[1,:]) #share of credit constrained agents (over all types)
     XX = [CC Sp]
     XX= XX[sortperm(XX[:, 1],rev=true), :]
     ca=cumsum(XX[:,2])

     nc= findall(x->x>Sh, ca)
     lim = nc[1] #number of credit constrained bins
     Th = sort(CC,rev=true)[lim] # the item nc of CC in decreasing order

     Id = Matrix(I,Nbin,Nbin)
     indcc=findall(x->x>=Th,CC)       #index of bin where fraction of CC is above Th. Those will be the bins where the agents are credit constrained
     indnc=findall(x->x<Th,CC)        #index of bin where fraction of CC is above Th. Those will be the bins where the agents are not credit constrained

     PP = copy(Id) #this is a diagonal matrix having 1 on the diagonal only if the history is not credit constrained
     for ii=1:length(indcc)
         PP[indcc[ii],indcc[ii]]= 0
     end
     PPc = Id - PP #this is the definition of Pc in the paper


     #############Constructing the xsis
     li = ls./Sp #per capita
     #li = ls
     Onev = ones(Float64,Nbin)

     ul_bin = fu(X)
     uc_bin = fup(X)
     ucc_bin = fupp(X)
     Mat2 = Id - β*R*(Transp')

     xsiu0 = ul./ul_bin
     xsiu1 = uc./uc_bin
     xsiu2 = ucc./ucc_bin

     xsiy = (lhsv./Sp)./((states[ytypeb].*li).^(1-τ))
     #xsiy = (lhsv)./((states[ytypeb].*li).^(1-τ))

     xsiuE = (inv(Mat2)*(resp))./(uc_bin)   #not per capita
     xsiv0 = (v0b) ./ fv.(li)

     xsiUc = (inv(Mat2)*(resp))./(uc_bin)   # per capita
     Res = (xsiUc.*uc_bin) -  β*R*(Transp')*(xsiUc.*uc_bin) - resp

     xsiv1 = (1-τ)*w.*xsiy.*(states[ytypeb].^(1-τ)).*(li.^(-τ)).*xsiu1.*uc_bin./(fvp.(li))  #xsil on labor supply equation with xsi1
     xsil = xsiv1

     #checks
     CC = zeros(6)
     CC[1] = sum(conso) - sum(dist.*gc) #conso
     CC[2] = sum(at) - sum(dist.*ga) #saving
     CC[3] = sum(at) - sum(ae2) #saving beginning and end
     CC[4] = maximum(Res) #euler equation on the bin
     CC[5] = maximum(xsil.*(li.^(1/ϕ)/χ) - (1-τ)*w.*xsiy.*(states[ytypeb].^(1-τ)).*(li.^(-τ)).*xsiu1.*uc_bin )
     CC[6] = maximum(at - Transp*ae)

     xsis = Xsis(xsiu0,xsiu1,xsiuE,xsiv0,xsiv1,xsiy,xsil)


     Proj = Projection(
     N,
     Nbin,
     R,
     w,
     A,
     B,
     states,
     Sp,
     at,
     ae,
     ls,
     conso,
     ul_bin,
     Transp,
     uc_bin,
     ucc_bin,
     Res,
     resp,
     xsis,
     Ntoi,
     CC,
     indnc,
     indcc,
     ytypeb,
     Li,
     Ri,
     Uc,
     Vl)
     return Proj,CC
 end
