#########################################################
 # Below the functions to solve Aiyagari through EGM
#########################################################


function interpEGM(
    pol::AbstractArray,
    grid::AbstractArray,
    x::T,
    na::Integer) where{T <: Real}
    np = searchsortedlast(pol,x)

    #if you give x (current asset) (out of the policy pol), gives ap (saving decision) out of the grid (grid)
    (np > 0 && np < na) ? np = np : #adjust indices if assets fall out of bounds
        (np == na) ? np = na-1 :
            np = 1

    ap_l,ap_h = pol[np],pol[np+1]
    a_l,a_h = grid[np], grid[np+1]
    ap = a_l + (a_h-a_l)/(ap_h-ap_l)*(x-ap_l)

    above =  ap > 0.0
    return above*ap,np #return either 0 or ap
end


function get_cEGM(
    pol::AbstractArray,
    Aggs::AggVarsEGM,
    CurrentAssets::AbstractArray,
    AiyagariModel::AiyagariModelEGM)

    @unpack params,aGrid,na,ns,states = AiyagariModel
    @unpack  γ,χ,ϕ,τ = params

    cpol = copy(pol)
    pol = reshape(pol,na,ns)
    for si = 1:ns
        for ai = 1:na
            asi = (si - 1)*na + ai
            F(c) = c-(Aggs.R*CurrentAssets[asi] +  Aggs.w*(states[si]* (χ*(1-τ)*Aggs.w *(states[si]^(1-τ))*(c^-γ))^(1/(1/ϕ+τ)))^(1-τ) - interpEGM(pol[:,si],aGrid,CurrentAssets[asi],na)[1] + Aggs.TT)
            csol = find_zero(F,1)
            cpol[asi] = csol
        end
    end
    return cpol
end


function EulerBackEGM(
    pol::AbstractArray,
    Aggs::AggVarsEGM,
    AiyagariModel::AiyagariModelEGM)

    @unpack params,na,ns,aGrid,Transv,states = AiyagariModel
    @unpack γ,β,χ,ϕ,τ = params

    aGridl = repeat(aGrid,ns)

    Trans = reshape(Transv,ns,ns)
    R,w,TT = Aggs.R,Aggs.w,Aggs.TT

    cp = get_cEGM(pol,Aggs,aGridl,AiyagariModel) #vector ns*na

    upcp = cp.^-γ
    Eupcp = copy(upcp)
    for ai = 1:na
        for si = 1:ns
            asi = (si-1)*na + ai
            Eupcp_sp = 0.0
            for spi = 1:ns
                aspi = (spi-1)*na + ai
                Eupcp_sp += Trans[spi,si]*upcp[aspi]
            end
            Eupcp[asi] = Eupcp_sp
        end
    end

    upc = β*Eupcp*R
    c = copy(upc)
    l = copy(upc)
    c = upc.^(-1.0/γ)
    apol = copy(upc)
    for ai = 1:na
        for si = 1:ns
            asi = (si-1)*na+ai
            l[asi]  = ((1-τ)*χ* w*(states[si]^(1-τ))*upc[asi])^(1/(1/ϕ+τ))
            apol[asi] = (aGridl[asi] + c[asi] - TT - w*(states[si]*l[asi])^(1-τ))/R
        end
    end
    #apol is a(a')
    #c is c(a')
    return apol
end


function SolveEGM(
    pol::AbstractArray,
    Aggs::AggVarsEGM,
    AiyagariModel::AiyagariModelEGM,
    tol = 1e-8)
    @unpack ns,na = AiyagariModel

    a = copy(pol)
    c = copy(pol)
    upc = copy(pol)
    l = copy(pol)

    for i = 1:10000
        a = EulerBackEGM(pol,Aggs,AiyagariModel)
        if (i-1) % 50 == 0
            test = abs.(a - pol)./(abs.(a) + abs.(pol))
            println("iteration: ",i," ",maximum(test))
            if maximum(test) < tol
                println("Solved in ",i," ","iterations")
                break
            end
        end
        pol = copy(a)
    end
    return pol
end


#########################################################
 # Below the functions  for the stationary distribution
#########################################################


function MakeTransMatEGM(pol,AiyagariModel,tmat)
    @unpack ns,na,aGrid,Transv = AiyagariModel
    pol = reshape(pol,na,ns)

    Trans = reshape(Transv,ns,ns)
    tmat = copy(tmat)
    tmat  .=0.0

    for a_i = 1:na
        for j = 1:ns
            x,i = interpEGM(pol[:,j],aGrid,aGrid[a_i],na)
            p = (aGrid[a_i] - pol[i,j])/(pol[i+1,j] - pol[i,j])
            p = min(max(p,0.0),1.0)
            sj = (j-1)*na
            for k = 1:ns
                sk = (k-1)*na
                tmat[sk+i+1,sj+a_i] = p * Trans[k,j]
                tmat[sk+i,sj+a_i] = (1.0-p) * Trans[k,j]
            end
        end
    end
    return tmat
end


function StationaryDistributionEGM(Mat,AiyagariModel)
    @unpack params,ns,na = AiyagariModel
    x = zeros(Float64,(na*ns,1))
    λ, x = powm!(Mat, ones(ns*na)/(ns*na), maxiter = 100000,tol = 1e-8)
    x = x/sum(x)
    return x
end


function steady(
    initialpol::AbstractArray,
    AiyagariModel::AiyagariModelEGM,
    R::T,
    w::T,
    TT::T,
    tol = 1e-8,maxn = 50)where{T <: Real}
    @unpack params,aGrid,ns,na,yvect = AiyagariModel

    tmat = zeros(eltype(initialpol),(na*ns,na*ns))
    pol = copy(initialpol)

    Aggs = AggVarsEGM(R,w,TT)
    pol = SolveEGM(pol,Aggs,AiyagariModel,tol)
    trans = MakeTransMatEGM(pol,AiyagariModel,tmat)
    D = StationaryDistributionEGM(trans,AiyagariModel)

    resE = copy(pol)
    resC = copy(pol)
    ga = copy(pol)  #policy rules as a function of beginning of period asset
    gl = copy(pol)  #policy rules as a function of beginning of period asset
    gc = copy(pol)  #policy rules as a function of beginning of period asset
    pol = reshape(pol,na,ns)
    for ai = 1:na
        for si = 1:ns
            x,i = interpEGM(pol[:,si],aGrid,aGrid[ai],na)
            asi = (si - 1)*na + ai
            F(c) = c-(R*aGrid[ai] + TT + w*(states[si]* (χ*(1-τ) *w*(states[si]^(1-τ))*c^-γ)^(1/(1/ϕ+τ)))^(1-τ) - x)
            gc[asi] = find_zero(F, 1)
            gl[asi]  = ((1-τ)*χ*w*(states[si]^(1-τ))*gc[asi]^-γ)^(1/(1/ϕ +τ))
            ga[asi] = x
            resC[asi] = gc[asi]-(R*aGrid[ai] + TT+ w*(states[si]* gl[asi])^(1-τ) - ga[asi])
        end
    end

    Rt       = 1/β
    A        = sum(D.*ga)
    Ltot     = (D'*(yvect.*gl))[1]
    K        = Ltot*( (Rt - (1-δ))/α )^(1/(α-1))
    resE     = gc.^-γ - β*R*trans'*(gc.^-γ)  #checking Euler equation

    Solution = AiyagariSolution(ga,gl,gc,D,trans,resE,R,w,TT,A,K,Ltot)
    return Solution
end
