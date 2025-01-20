########################################################################################
####### PACKAGES TO RUN THE CODE
########################################################################################

using LinearAlgebra
using Parameters
using IterativeSolvers
using FastGaussQuadrature
using ForwardDiff
using QuantEcon
using Plots
using Arpack
using BenchmarkTools
using JLD2
using MATLAB
using Plots; pyplot();
using SparseArrays;
using MAT
using Statistics
using Roots
using Interpolations
using NLsolve
using CSV
using DataFrames


########################################################################################
####### STRUCUTURES FOR THE CODE
########################################################################################

struct AiyagariParametersEGM{T <: Real}
    β::T
    α::T
    δ::T
    γ::T
    χ::T
    ϕ::T
    κ::T #labor tax level
    τ::T #labor tax progressivity
    tc::T #consumption tax
    tk::T #capital tax
    TT::T #transfer
end


struct Aiyagarifunction
    fu::Function
    fuinv::Function
    fup::Function
    fupp::Function
    fv::Function
    fvp::Function
    fvpp::Function
end


struct AiyagariModelEGM{T <: Real,I <: Integer}
    params::AiyagariParametersEGM{T}
    Functions::Aiyagarifunction
    aGrid::Array{T,1} #policy grid
    na::I #number of grid points in policy function
    states::Array{T,1} #earning states
    ns::I #number of states
    Transv::Array{T,1}
    disty::Array{T,1} #distribution y size is ns
    yvect::Array{T,1} #vector of y size is ns*na
end


struct AiyagariSolution{T <: Real}
    ga::Array{T,1} #policy grid
    gl::Array{T,1}
    gc::Array{T,1} #number of grid points in policy function
    dist::Array{T}
    Mat::Matrix{T}
    resE::Array{T}
    R::T
    w::T
    TT::T
    A::T
    K ::T
    Ltot :: T
end


mutable struct AggVarsEGM{S <: Real,T <: Real}
    R::S
    w::T
    TT::T
end


struct Xsis{T <: Real}
    xsiu0::Array{T,1} #correcting paremeter for utility in level
    xsiu1::Array{T,1} #correcting paremeter for u'
    xsiuE::Array{T,1} #correcting paremeter for Euler equation
    xsiv0::Array{T,1} #correcting paremeter for disutility of labor
    xsiv1::Array{T,1} #correcting paremeter for v'
    xsiy::Array{T,1} #correcting paremeter for labor
    xsil::Array{T,1} #correcting paremeter for v'
end


## Be careful the following structure store information by bin (and not per capita). Divide by the size of the bin to have per capita terms
struct Projection{T <: Real,I <: Integer}
    N::I
    Nbin::I #number of bins
    R::T
    w::T
    A::T
    B::T
    states::Array{T,1}
    Sp::Array{T,1}  #size of each bin
    abp::Array{T,1} #not per capita
    aep::Array{T,1} #not per capita
    lp::Array{T,1}  #not per capita
    cp::Array{T,1}  #not per capita
    ul::Array{T,1}  #vector of utilities by bin
    Matab::SparseMatrixCSC{T,I}
    uc::Array{T,1} #vector of marginal utilities by bin
    ucc::Array{T,1} #vector of derivative of marginal utility by bin
    Res::Array{T,1} #check that resb is small
    resp::Array{T,1} #lagrange multiplier on credit constraint by bin (not per capita)
    xsis:: Xsis
    Ntoi::Array{I,1} #for euler equations
    CC::Array{T,1}
    indnc::Array{I,1} #index of credit constrained histories
    indcc::Array{I,1} #Index of credit constrained histories
    ytype::Array{I,1}
    Li::Array{T,1}
    Ri::Array{T,1}
    Uc::Array{T,1}
    Vl::Array{T,1}
end


## Variables of the planner
struct Planner{T <: Real}
    lambdac::Array{T,1}
    lambdal::Array{T,1}
    lambdact::Array{T,1}
    psi::Array{T,1}
    weight::Array{T,1} #weights for all bins
    weightp::Array{T,1} #weights for productivity only
    weightb::Array{T,1} #weights for all the bins to be used in MATLAB
    mu::T
end


## Function to compute the Gini
function gini2(wagedistarray)
    Swages = cumsum(wagedistarray[:,1].*wagedistarray[:,2])
    Gwages = Swages[1]*wagedistarray[1,2] +
                sum(wagedistarray[2:end,2] .*
                        (Swages[2:end]+Swages[1:end-1]))
    return 1 - Gwages/Swages[end]
end


function wealthD(ww::AbstractArray,NN::Integer)
ww2 = sortslices(ww, dims=1)
fracw = cumsum(ww2[:,1].*ww2[:,2])
fracw =fracw ./fracw[end]
fracS = cumsum(ww2[:,2])
plot(fracS,fracw)
vect = zeros(Float64,NN,2)
is  = 0
for i=1:NN
    is = maximum(findall(x->x<=i/NN,fracS))
    vect[i,1] = fracS[is] #- vect[i-1,1]
    vect[i,2] = fracw[is] #- vect[i-1,2]
end
vecf = zeros(Float64,NN,2)
vecf[1,:] .= vect[1,:]
for i=2:NN
    vecf[i,1] = vect[i,1] - vect[i-1,1]
    vecf[i,2] = vect[i,2] - vect[i-1,2]
end
return  vecf
end


## Exponential grid
function grid_fun(a_min,a_max,na, pexp)
    x = range(a_min,step=0.5,length=na)
    grid = a_min .+ (a_max-a_min)*(x.^pexp/maximum(x.^pexp))
    return grid
end


########################################################################################
####### RUNNING THE MAIN CODE TO OBTAIN THE STEADY STATE, REPRODUCE THE TAX SYSTEM, AND 
####### OBTAIN THE INCOME AND WEALTH DISTRIBUTION
########################################################################################

include("parameters_France_final.jl") #file with the main parameters to obtain the steady state
include("Aiyagari_solve_endo.jl") #function to solve the heterogenous agent model in steady state

Solution  = steady(pol0,Modeling,R,w,TT,1e-8)


## Checking solution
@unpack A,K,Ltot,gc,gl,dist,ga = Solution
@unpack yvect = Modeling
Ltotb = (dist'*((yvect.*gl).^(1-τ)))[1] #aggregate for the labor in the budget constraint of consumption
K = Ltot*( (Rt - (1-δ))/α )^(1/(α-1))
B  = (1+tc)*A - K
Y = K^α*Ltot^(1-α)
Ctot = w*Ltotb + (R-1)*(A) + TT # = sum(gc.*dist)
G = Y- δ*K - Ctot
UG = 0
#UG = (1/(χG))*(G^(1-θ) -1)./(1-θ)
#UGp = (1/(χG))*(G)^(-θ)
UGp = 0


## Various additional checks
G2 = Y + δ*(B - tc*A) - (R-1 +δ)*(A) - w*Ltotb - TT
CC1 = zeros(2)
CC1[1] = (Ctot - sum(gc.*dist))/Ctot
CC1[2] = G - G2


## Distribution
@unpack aGrid = Modeling
aGridl = repeat(aGrid,ns)
D2 = copy(dist)
D2[1] = dist[1]
sum(D2)
ww = [aGridl D2]
ww2 = sortslices(ww, dims=1)
Gini = gini2(ww2)

ww = [gc D2] #labor supply I don't consider wage, as does change Gini
ww2 = sortslices(ww, dims=1)
Gini_c = gini2(ww2)

#vW = wealthD(ww2,5)


## Credit Constrained
distv = reshape(dist,na,ns) #distribution of wealth
Sh = sum(distv[1,:]) #share of credit constrained agents (over all types)
format(x) = round(x*10^2)/10^2


## Calculating the income after tax and transfers
income_at =Array{Float64,1}(undef, ns*na)
for si = 1:ns
    for ai = 1:na
        asi = (si - 1)*na + ai
         income_at[asi] = R*ga[asi] -aGridl[asi] + w*(states[si]*gl[asi])^(1-τ) + TT
    end
end


## Calculating the income before taxes and transfers
income_bt =Array{Float64,1}(undef, ns*na)
for si = 1:ns
    for ai = 1:na
        asi = (si - 1)*na + ai
         income_bt[asi] = (Rt-1)*aGridl[asi]   + wt*states[si]*gl[asi]
    end
end


## Gini coefficient for Income After Tax
wwat = [income_at D2]
ww2at = sortslices(wwat, dims=1) #put in order the matrix according to the row, which represents the income in order to calculate the Gini index
Gini_income_at = gini2(ww2at)
#vWat = wealthD(ww2at,5)


## Gini coefficient for Income Before Tax and Transfers
wwbt = [income_bt D2]
ww2bt = sortslices(wwbt, dims=1) #put in order the matrix according to the row, which represents the income in order to calculate the Gini index
Gini_income_bt = gini2(ww2bt)
#vWbt = wealthD(ww2bt,5)


## Summary statistics for taxes
rtk = (tk*(Rt-1)*(A))/Y #total capital tax as % of GDP
rtl =(wt*Ltot - (κ*wt^(1-τ)).*Ltotb)/Y  #total labor tax as % of GDP
rtc=(tc*Ctot)/Y  #total consumption tax as % of GDP
TPO = rtk + rtl + rtc #total tax as % of GDP

cpol = reshape(gc,na,ns)
grida = reshape(aGridl,na,ns)
dist = reshape(D2,na,ns)
dc = cpol[2:na,:] .- cpol[1:na-1,:]
da = grida[2:na,:].- grida[1:na-1,:]
mpc = dc./da
ampc = (sum(mpc.*(dist[1:na-1,:])))


## Results of the Steady State Macroeconomic allocations and Inequality Results
println("B/Y: ",format(B/(Y*4)), " Gini: ", format(Gini), " G/Y: ",format(G/Y), " C/Y: ",format(Ctot/Y), " I/Y: ",format(δ*K/Y), " Tr/Y: ",format(TT/Y)," Ltot: ",format(Ltot), " MPC: ", format(ampc), " TPO: ", format(TPO), " Share : ",format(Sh))
CC1


########################################################################################
####### TRUNCATION METHOD
########################################################################################

include("Projection_truncation_endo.jl") #code to use the Truncation Method
N  = 5 #Truncation length
proj,CC = Project_plan(N,Modeling,Solution)
println(CC);


## Computing the Weights
positive = true
include("Function_pareto_weights_G_optimal.jl")
planner,CC,J1,J2,J3,J4,J5,M7 = Weight(proj,Modeling,Solution)
planner.weight
plot(planner.weightp,line = (3,:solid, [:blue]), label = "Weights")
plot!(planner.weight, label = "Weights Marginal Utility ")
plot!(10*disty,legend=:topleft, line = (2, :dash, [:darkgreen]),fill = (0, :lightgreen),fillalpha = 0.3, label = "Distribution")

z = log.(states)
function f!(F, x)
    F[1] = sum(J2*exp.(x[1]*z.^0+x[2]*z.^1+x[3]*z.^2))
    F[2] = sum(J3*exp.(x[1]*z.^0+x[2]*z.^1+x[3]*z.^2))
    F[3] = sum(J5*exp.(x[1]*z.^0+x[2]*z.^1+x[3]*z.^2)) .-1
end


## Finding the coefficients for the Parametric Weights
sol = nlsolve(f!, [0.; 0.; 0.;],autodiff = :forward);
a = sol.zero[1]
b = sol.zero[2]
c = sol.zero[3]
weightparam =exp.(a*z.^0+b*z.^1 + c*z.^2)
weightparam = weightparam/sum(disty.*weightparam)

weights = planner.weightp./mean(planner.weightp)
weights_marginal = planner.weight./mean(planner.weight)
weights_param = weightparam./mean(weightparam)


## Plotting the Non-parametric Weights and the Marginal Utility Weights
plot(weights,line = (3,:solid, [:blue]), label = "Weights")
plot!(weights_marginal, label = "Weights Marginal Utility ")
plot!(10*disty,legend=:topleft, line = (2, :dash, [:darkgreen]),fill = (0, :lightgreen),fillalpha = 0.3,dpi = 300, label = "Distribution")
ylims!(0,4)
savefig("Weights_france.png")


## Plotting the Non-parametric Weights and the Parametric ones
plot(weights,line = (3,:solid, [:blue]), label = "Weights")
#plot!(weights_marginal, label = "Weights Marginal Utility ")
plot!(weights_param, line = (2, :dash, [:black]), label = "Weights Parametric ")
plot!(10*disty,legend=:topleft, line = (2, :dash, [:darkgreen]),fill = (0, :lightgreen),fillalpha = 0.3,dpi = 300, label = "Distribution")
ylims!(0,5.5)
savefig("Weights_total_france.png")


## Plotting the Parametric Weights
plot(weights_param, line = (3,:solid, [:blue]), label = " SWF Weights ",legend_frame=false)
plot!(10*disty,legend=:topleft,legend_frame=false,line = (2, :dash, [:darkgreen]),fill = (0, :lightgreen),fillalpha = 0.3,dpi = 300, label = "Distribution")
ylims!(0,5.5)
savefig("Weights_param_France.png")


########################################################################################
####### COMPUTING GSMWW
########################################################################################

x = cumsum(disty*100)
uc = 1.0./weights_marginal #weights_marginal are constructed as 1/uc
y = weights_param.*uc
yn = y/y[5]


## For IOWD we consider the average weight per decile
ycum = cumsum(disty.*y)
itp_g = interpolate((x,), ycum, Gridded(Linear()))
ext_itp_g = extrapolate(itp_g, Line())
IOWd =zeros(Float64,5)
for i=1:5
    IOWd[i]=ext_itp_g(i*20)-ext_itp_g((i-1)*20+2)
end
IOWd=IOWd/mean(IOWd)


## For IOW we plot the weight for each productivity level
itp_w = interpolate((x,), yn, Gridded(Linear()))
ext_itp_w = extrapolate(itp_w, Line())
IOW =zeros(Float64,10)
for i=1:10
    IOW[i] = ext_itp_w(i*10)
end
plot(x[1:end],IOW[1:end])


## Plotting GSMWW
plot(IOWd,legend=:false,line = (3,:solid, [:blue]))
savefig("GSMWW_France.png")


## Computing the Matrix to decentralize
df = CSV.read("data_france.csv", DataFrame)
x = Float64.(df[:,1])
y = df[:,2]
Wage2007 = 31093 # in2007 euros, average family income
x_new = Modeling.states*Wage2007
x[end] = x_new[end]

itp = interpolate((x,), y, Gridded(Linear()))
#x[1] = x_new[1] # to be sure to be in the grid
ext_itp = extrapolate(itp, Line())


WeightP = zeros(Float64,length(x_new))
WeightP = [ext_itp(xi) for xi in x_new]
WeightP = WeightP*length(WeightP)/sum(WeightP) # to have it around 1


## Plotting the Participation Rate as a Fraction of the Average Wage
x_plot = 100*x_new/Wage2007
extra= [ext_itp(xi) for xi in x_new]
plot(x_plot,extra,legend=false,xlims=(30,300),line = (3,:solid, [:blue]))
savefig("participation_france.png")


## Normalization with political weights, for the weighted sum to be 1
WeightS = weights_param/sum(weights_param.*Modeling.disty.*WeightP)


## Computing the matrices
include("Function Weight_Matrix.jl")
WM,WM2,WM2s = Weight_Matrix(Modeling,WeightP,WeightS)


## Plotting the Individual Welfare Functions (IWFs)
plot([1:length(WeightS)],[1:length(WeightS)],WM2',st=:surface,
xlabel= "Actual productivity",ylabel ="Considered productivty",zlabel="Weights",legend=false,alpha=0.7,camera=(-26, 17),zlims=(0,1.6))
savefig("Weights_france_surf.png")


plot([1:length(WeightS)],[1:length(WeightS)],WM2,st=:surface,
xlabel= "Actual productivity",ylabel ="Considered productivty",zlabel="Weights",legend=false,alpha=0.7,camera=(-20, 20))


## Plotting the Self-seeking Weights 
plot([1:length(WeightS)],diag(WM2),xlabel= "Actual productivity",ylabel ="Self-seeking Weights",legend=false,alpha=0.7)
savefig("Weights_self_france.png")


## Summary statistics of the Non-Parametric Weights
quantile(weights, [0.1,0.5, 0.9])
round(mean(weights),sigdigits=3)
round(std(weights),sigdigits=3)
round(minimum(weights),sigdigits=3)
round(maximum(weights),sigdigits=3)


########################################################################################
####### STUDYING THE IMPACT OF DIFFERENT PARAMETERS, TAX SYSTEMS, AND 
####### INCOME PROCESS 
########################################################################################

M8= zeros(proj.Nbin,ns)
for i=1:proj.Nbin
    M8[i,proj.ytype[i]]=1
end

Utility = (M8'*(proj.Sp.*proj.Uc)./disty) - (M8'*(proj.Sp.*proj.Vl)./disty)
Labor_Income = (M8'*(proj.Sp.*proj.Li)./disty)
Capital_Income = (M8'*(proj.Sp.*proj.Ri)./disty)


## Government spending as a % of GDP
Gp = (G./Y)


## Define a function to run the computations for different parameter files and save results
function run_computation_and_save(parameters_file, suffix)
    
    global κ  
    global τ 
    global aGrid
    global yvect
    include(parameters_file)

    #Gp = (G ./ Y)
    κ1 = κ:κ
    τ1 = τ:τ
    err_κ = 0.01
    err_τ = 0.01
    tol_G = 0.001

    R1 = zeros(Float64, length(κ1) * length(τ1))
    w1 = zeros(Float64, length(κ1) * length(τ1))

    for i = 1:length(κ1)
        for j = 1:length(τ1)
            asi = (i - 1) * length(τ1) + j
            R1[asi] = (1 + (1 - tk) * (Rt - 1)) * (1 + tc) / (1 + tc)
            w1[asi] = (κ1[i] * wt^(1 - τ1[j])) / (1 + tc)
        end
    end

    κ2 = kron(κ1[1:end], ones(length(τ1), 1))
    τ2 = kron(ones(length(κ1), 1), τ1[1:end])

    A1 = zeros(Float64, length(κ2))
    K1 = zeros(Float64, length(κ2))
    Ltot1 = zeros(Float64, length(κ2))
    Ltotb1 = zeros(Float64, length(κ2))
    B1 = zeros(Float64, length(κ2))
    Ctot1 = zeros(Float64, length(κ2))
    Y1 = zeros(Float64, length(κ2))
    G1 = zeros(Float64, length(κ2))
    G21 = zeros(Float64, length(κ2))
    GY = zeros(Float64, length(κ2))

    gc1 = zeros(Float64, na * ns, length(κ2))
    gl1 = zeros(Float64, na * ns, length(κ2))
    dist1 = zeros(Float64, na * ns, length(κ2))

    dist_G = ones(Float64, length(κ2))
    tolG = 0.005 * ones(Float64, length(κ2))

    for k = 1:length(dist_G)
        κ = κ2[k]
        τ = τ2[k]
        params = AiyagariParametersEGM(β, α, δ, γ, χ, ϕ, κ, τ, tc, tk, TT)
        pol0 = repeat(aGrid, ns)
        Modeling = AiyagariModelEGM(params, Functions, aGrid, na, states, ns, Transv, disty, yvect)

        Solution = steady(pol0, Modeling, R1[k], w1[k], TT, 1e-8)
        @unpack A, K, Ltot, gc, gl, dist = Solution
        A1[k], K1[k], Ltot1[k], gc1[:, k], gl1[:, k], dist1[:, k] = A, K, Ltot, gc, gl, dist

        Ltotb1[k] = (dist1[:, k]' * ((yvect .* gl1[:, k]) .^ (1 - τ)))[1]  #aggregate labor for the budget of the economy 
        K1[k] = Ltot1[k] * ((Rt - (1 - δ)) / α)^(1 / (α - 1))
        B1[k] = (1 + tc) * A1[k] - K1[k]
        Y1[k] = K1[k]^α * Ltot1[k]^(1 - α)
        Ctot1[k] = w1[k] * Ltotb1[k] + (R1[k] - 1) * (A1[k]) + TT # = sum(gc1[:, k] .* dist1[:, k])
        G1[k] = Y1[k] - δ * K1[k] - Ctot1[k]
        G21[k] = Y1[k] + δ * (B1[k] - tc * A1[k]) - (R1[k] - 1 + δ) * (A1[k]) - w1[k] * Ltotb1[k] - TT

        GY[k] = G1[k] / Y1[k]
        dist_G[k] = abs(GY[k] - Gp)
    end

    distG, ind = findmin(dist_G)
    κ = κ2[ind]
    τ = τ2[ind]
    params = AiyagariParametersEGM(β, α, δ, γ, χ, ϕ, κ, τ, tc, tk, TT)
    pol0 = repeat(aGrid, ns)
    Modeling = AiyagariModelEGM(params, Functions, aGrid, na, states, ns, Transv, disty, yvect)

    Solution = steady(pol0, Modeling, R1[ind], w1[ind], TT, 1e-8)

    include("Projection_truncation_endo.jl")
    N = 5 # truncation length
    proj, CC_proj = Project_plan(N, Modeling, Solution)
    println(CC_proj)

    include("Function_pareto_weights_G_optimal.jl")
    planner, CC_planner, J1, J2, J3, J4, J5 = Weight(proj, Modeling, Solution)
    weights = planner.weight

    z = log.(states)
    function f!(F, x)
        F[1] = sum(J2 * exp.(x[1] * z.^0 + x[2] * z.^1 + x[3] * z.^2))
        F[2] = sum(J3 * exp.(x[1] * z.^0 + x[2] * z.^1 + x[3] * z.^2))
        F[3] = sum(J5 * exp.(x[1] * z.^0 + x[2] * z.^1 + x[3] * z.^2)) - 1
    end

    sol = nlsolve(f!, [0.; 0.; 0.], autodiff = :forward)
    a, b, c = sol.zero

    weightparam = exp.(a * z.^0 + b * z.^1 + c * z.^2)
    weightparam /= sum(disty .* weightparam)

    weights_p = planner.weightp ./ mean(planner.weightp)
    weights_marginal = planner.weight ./ mean(planner.weight)
    weights_param = weightparam ./ mean(weightparam)

    M8 = zeros(proj.Nbin, ns)
    for i = 1:proj.Nbin
        M8[i, proj.ytype[i]] = 1
    end

    Utility = (M8' * (proj.Sp .* proj.Uc) ./ disty) - (M8' * (proj.Sp .* proj.Vl) ./ disty)
    Labor_Income = (M8' * (proj.Sp .* proj.Li) ./ disty)
    Capital_Income = (M8' * (proj.Sp .* proj.Ri) ./ disty)

    κo = format(κ2[ind])

    weights = [weights, weights_p, weights_marginal, weights_param, disty]

    save("weights_$suffix.jld", "weights", weights)
end


## Run the function for each parameter file and save results
run_computation_and_save("parameters_France_final_beta.jl", "beta") #French parameters with US beta
run_computation_and_save("parameters_France_final_beta_taxes.jl", "beta_taxes") #French parameters with US beta and taxes
run_computation_and_save("parameters_France_final_beta_inc_taxes.jl", "beta_inc_taxes") #French parameters with US beta, taxes, and income process


## Loading the results for plotting
using Plots
weights_beta = load("weights_beta.jld", "weights")[2]
weights_beta_taxes = load("weights_beta_taxes.jld", "weights")[2]
weights_beta_inc_taxes = load("weights_beta_inc_taxes.jld", "weights")[2]


## Plotting the results for the Non-Parametric Weights
plot(weights_beta, line = (3, :dot, [:red]), label = "Weights with US β")
plot!(weights_beta_taxes, line = (3, :dash, [:lightred]), label = "Weights with US β + Taxes")
plot!(weights_beta_inc_taxes, line = (3, :solid, [:blue]), label = "Weights with US β + Taxes + Income Process")
plot!(10 * disty, legend = :topleft, line = (2, :dash, [:darkgreen]), fill = (0, :lightgreen), fillalpha = 0.3, dpi = 300, label = "Distribution")
ylims!(0, 5)
savefig("Weights_comparison_total_France_2.png")


## Loading the results for plotting
weights_beta = load("weights_beta.jld", "weights")[4]
weights_beta_taxes = load("weights_beta_taxes.jld", "weights")[4]
weights_beta_inc_taxes = load("weights_beta_inc_taxes.jld", "weights")[4]


## Plotting the results for the Parametric Weights
plot(weights_beta, line = (3, :dot, [:red]), label = "Weights with US β")
plot!(weights_beta_taxes, line = (3, :dash, [:lightred]), label = "Weights with US β + Taxes")
plot!(weights_beta_inc_taxes, line = (3, :solid, [:blue]), label = "Weights with US β + Taxes + Income Process")
plot!(10 * disty, legend = :topleft, line = (2, :dash, [:darkgreen]), fill = (0, :lightgreen), fillalpha = 0.3, dpi = 300, label = "Distribution")
ylims!(0, 5)
savefig("Weights_comparison_total_France_param_2.png")


########################################################################################
####### STUDYING THE IMPACT OF CHANGING THE FISCAL SYSTEM
########################################################################################

M8= zeros(proj.Nbin,ns)
for i=1:proj.Nbin
    M8[i,proj.ytype[i]]=1
end

Utility = (M8'*(proj.Sp.*proj.Uc)./disty) - (M8'*(proj.Sp.*proj.Vl)./disty)
Labor_Income = (M8'*(proj.Sp.*proj.Li)./disty)
Capital_Income = (M8'*(proj.Sp.*proj.Ri)./disty)


## Government spending as a % of GDP
Gp = (G./Y)


## Using the Parameters of US
include("parameters_US_final.jl")
κ1 = 0.64:0.001:0.65
τ1 = τ:τ
#τ1 = 0.78:0.001:0.79
err_κ = 0.01
err_τ = 0.01
tol_G = 0.001


R1 = zeros(Float64,length(κ1)*length(τ1))
w1 = zeros(Float64,length(κ1)*length(τ1))


for i =1:length(κ1)
    for j = 1:length(τ1)
        asi = (i - 1)*length(τ1) + j
        R1[asi] = (1 + (1 - tk)*(Rt-1))*(1+tc)./(1+tc)
        w1[asi] = (κ1[i]*wt^(1 .-τ1[j]))./(1+tc)
    end
end

κ2 = kron(κ1[1:end],ones(length(τ1),1))
τ2 = kron(ones(length(κ1),1), τ1[1:end])

A1 = zeros(Float64,length(κ2))
K1 = zeros(Float64,length(κ2))
Ltot1 = zeros(Float64,length(κ2))
Ltotb1 = zeros(Float64,length(κ2))
B1 = zeros(Float64,length(κ2))
Ctot1 = zeros(Float64,length(κ2))
Y1 = zeros(Float64,length(κ2))
G1 = zeros(Float64,length(κ2))
G21 = zeros(Float64,length(κ2))
GY = zeros(Float64,length(κ2))

gc1 = zeros(Float64,na*ns,length(κ2))
gl1 = zeros(Float64,na*ns, length(κ2))
dist1 = zeros(Float64,na*ns, length(κ2))

dist_G = ones(Float64,length(κ2))
tolG = 0.005*ones(Float64,length(κ2))


for k = 1:length(dist_G)

    #dist_G = ones(Float64,length(tk))
    #tolG = 0.005*ones(Float64,length(tk))

    #while dist_G[k] > tolG[k]
    κ = κ2[k]
    τ = τ2[k]
    params = AiyagariParametersEGM(β,α,δ,γ,χ,ϕ,κ,τ,tc,tk,TT)
    pol0   = repeat(aGrid,ns)
    Modeling  = AiyagariModelEGM(params, Functions, aGrid,na,states,ns,Transv,disty,yvect)

    Solution = steady(pol0,Modeling,R1[k],w1[k],TT,1e-8)
    @unpack A,K,Ltot,gc,gl,dist = Solution
    A1[k],K1[k],Ltot1[k],gc1[:,k],gl1[:,k],dist1[:,k] = A,K,Ltot,gc,gl,dist

    Ltotb1[k] = (dist1[:,k]'*((yvect.*gl1[:,k]).^(1-τ)))[1] # aggregate for the budget of
    K1[k] = Ltot1[k]*( (Rt - (1-δ))/α )^(1/(α-1))
    B1[k]  = (1+tc)*A1[k] - K1[k]
    Y1[k] = K1[k]^α*Ltot1[k]^(1-α)
    Ctot1[k] = w1[k]*Ltotb1[k] + (R1[k]-1)*(A1[k]) + TT # = sum(gc1[:,k].*dist1[:,k])
    G1[k] = Y1[k]- δ*K1[k] - Ctot1[k]
    ## various check
    G21[k] = Y1[k] + δ*(B1[k] - tc*A1[k]) - (R1[k]-1 +δ)*(A1[k]) - w1[k]*Ltotb1[k] - TT

    GY[k] = G1[k]/Y1[k]
    dist_G[k] =abs.(GY[k] .- Gp)
end


distG, ind = findmin(dist_G)
κ = κ2[ind]
τ=τ2[ind]
params = AiyagariParametersEGM(β,α,δ,γ,χ,ϕ,κ,τ,tc,tk,TT)
pol0   = repeat(aGrid,ns)
Modeling  = AiyagariModelEGM(params, Functions, aGrid,na,states,ns,Transv,disty,yvect)

Solution  = steady(pol0,Modeling,R1[ind],w1[ind],TT,1e-8)


## Truncation Method
include("Projection_truncation_endo.jl")
N  = 5 #Truncation length
proj_fr,CC_fr = Project_plan(N,Modeling,Solution)
println(CC_fr);


## Computing the weights
positive = true
include("Function_pareto_weights_G_optimal.jl")
planner_fr,CC_fr,J1_fr,J2_fr,J3_fr,J4_fr,J5_fr = Weight(proj_fr,Modeling,Solution)
planner_fr.weight

z = log.(states)
function f!(F, x)
    F[1] = sum(J2_fr*exp.(x[1]*z.^0+x[2]*z.^1+x[3]*z.^2))
    F[2] = sum(J3_fr*exp.(x[1]*z.^0+x[2]*z.^1+x[3]*z.^2))
    F[3] = sum(J5_fr*exp.(x[1]*z.^0+x[2]*z.^1+x[3]*z.^2)) .-1
end

sol = nlsolve(f!, [0.; 0.; 0.;],autodiff = :forward);
a = sol.zero[1]
b = sol.zero[2]
c = sol.zero[3]
weightparam_fr =exp.(a*z.^0+b*z.^1 + c*z.^2)
weightparam_fr = weightparam_fr/sum(disty.*weightparam_fr)


weights_fr = planner_fr.weightp./mean(planner_fr.weightp)
weights_marginal_fr = planner_fr.weight./mean(planner_fr.weight)
weights_param_fr = weightparam_fr./mean(weightparam_fr)


M8_fr = zeros(proj_fr.Nbin,ns)
for i=1:proj_fr.Nbin
    M8_fr[i,proj_fr.ytype[i]]=1
end


Utility_fr = (M8_fr'*(proj_fr.Sp.*proj_fr.Uc)./disty) - (M8_fr'*(proj_fr.Sp.*proj_fr.Vl)./disty)
Labor_Income_fr = (M8_fr'*(proj_fr.Sp.*proj_fr.Li)./disty)
Capital_Income_fr = (M8_fr'*(proj_fr.Sp.*proj_fr.Ri)./disty)


κo=format(κ2[ind])


## Plotting the Weights
plot(weights,line = (3,:solid, [:blue]), label = "Weights")
plot!(weights_fr,line = (3,:dash, [:red]), label = "Weights with US system")
plot!(weights_marginal, line = (2,:solid, [:lightblue]), label = "Weights Marginal Utility")
plot!(weights_marginal_fr, line = (2,:dash, [:lightred]), label = "Weights Marginal Utility with US system")
plot!(10*disty,legend=:topleft, line = (2, :dash, [:darkgreen]),fill = (0, :lightgreen),fillalpha = 0.3,dpi=300, label = "Distribution")
ylims!(-0.2,5.0)
#savefig("Weights_france_us.png")

plot(weights,line = (3,:solid, [:blue]), label = "Weights")
plot!(weights_fr,line = (3,:dash, [:red]), label = "Weights \$κ =\$ $κo")
plot!(Labor_Income, line = (2,:solid, [:lightblue]), label = "Labor Income ")
plot!(Labor_Income_fr, line = (2,:dash, [:lightred]), label = "Labor Income \$κ =\$ $κo")
plot!(10*disty,legend=:topright, line = (2, :dash, [:darkgreen]),fill = (0, :lightgreen),fillalpha = 0.3,dpi=300, label = "Distribution")

plot(weights,line = (3,:solid, [:blue]), label = "Weights")
plot!(weights_fr,line = (3,:dash, [:red]), label = "Weights \$κ =\$ $κo")
plot!(Capital_Income, line = (2,:solid, [:lightblue]), label = "Capital Income ")
plot!(Capital_Income_fr, line = (2,:dash, [:lightred]), label = "Capital Income \$κ =\$ $κo")
plot!(10*disty,legend=:topright, line = (2, :dash, [:darkgreen]),fill = (0, :lightgreen),fillalpha = 0.3,dpi=300, label = "Distribution")

l = @layout [a ; b c;d]
p1 = plot(weights_fr .- weights, line = (2,:solid, [:blue]), label = "Δ Weights", ylims = (-0.25,0.15), legend=:bottomright )
p2 = plot(Utility_fr .- Utility, line = (2,:solid, [:blue]), label = "Δ Utility", ylims =  (-0.25,0.15), legend=:bottomright)
p3 = plot(Labor_Income_fr .- Labor_Income, line = (2,:solid, [:blue]), label = "Δ Labor Income", ylims = (-0.25,0.15), legend=:bottomright)
p4 = plot(Capital_Income_fr .- Capital_Income, line = (2,:solid, [:blue]), label = "Δ Capital Income", ylims = (-0.25,0.15) , legend=:bottomright )
plot(p1, p2, p3,p4,dpi = 300, title=["Weights " "Utility" "Labor Inc. " "Capital Inc. "])
#savefig("Dif_Weights_france_us.png")


## Plotting the difference between the original calibration and the one using the US Tax System
l = @layout [a ; b c;d]
p1 = plot(weights_param_fr .- weights_param, line = (2,:solid, [:blue]), label = "Δ Weights", ylims = (-0.25,0.15), legend=:bottomright )
p2 = plot(Utility_fr .- Utility, line = (2,:solid, [:blue]), label = "Δ Utility", ylims = (-0.25,0.15), legend=:bottomright)
p3 = plot(Labor_Income_fr .- Labor_Income, line = (2,:solid, [:blue]), label = "Δ Labor Income", ylims = (-0.25,0.15), legend=:bottomright)
p4 = plot(Capital_Income_fr .- Capital_Income, line = (2,:solid, [:blue]), label = "Δ Capital Income", ylims = (-0.25,0.15), legend=:bottomright )
plot(p1, p2, p3,p4,dpi = 300, title=["Weights " "Utility " "Labor Inc. " "Capital Inc. "])
savefig("Dif_Weights_france_us_param.png")


########################################################################################
####### COMPUTING THE WEIGHTS BY HISTORY
########################################################################################

## Running the main code
include("parameters_France_final.jl")
include("Aiyagari_solve_endo.jl")

Solution  = steady(pol0,Modeling,R,w,TT,1e-8)


positive = true
include("Function_pareto_weights_Bin.jl") #function to obtain the results by bin
planner_bin,CC_bin,J1_bin,J2_bin,J3_bin,J4_bin,J5_bin = Weight(proj,Modeling,Solution)

weights_bin = planner_bin.weightp./mean(planner_bin.weightp)


## Normalizing the wages and the Wealth
wbar = 13.23 #average Annual Wages
ls = 8760*proj.lp./proj.Sp
wtn = (wbar*sum(proj.Sp.*ls))./sum(proj.Sp.*proj.xsis.xsiy.*(ls.*states[proj.ytype]).^(1-τ))
wage = wtn.*(ls.*states[proj.ytype]).^(1-τ)


sum(proj.Sp.*wage) #this gives the average wage


wealth = (proj.abp)./(proj.Sp)
GDPUS= 30332.328172
An = (K/(Y*4) + B/(Y*4))*GDPUS
unit = An./sum(proj.Sp.*wealth)
wealth = (proj.abp.*unit)./(proj.Sp)
wealth = wealth .-minimum(wealth)


using Plots; pyplot()
font = Plots.font("Times New Roman", 10)


plot(wealth./1000,states[proj.ytype],weights_bin,st=:surface, legend = false, zlabel=("Pareto weights"), zguidefontrotation=90, yguidefontrotation=-15, xguidefontrotation=45, dpi=300 ,camera=(-60,30))
zlims!(0.0,5.0)
xlabel!("Wealth in k of €")
ylabel!("Productivity level")
#savefig("weights_france_wealth.png")


plot(wealth./1000,wage./1000,weights_bin,st=:surface, legend = false, zlabel=("Pareto weights"), zguidefontrotation=90, yguidefontrotation=-15, xguidefontrotation=45, dpi=300 , camera=(-60,30))
zlims!(0.0,5.0)
xlabel!("Wealth in k of €")
ylabel!("Average wage per year in k of € ")
#savefig("weights_france_wealth_income.png")


wealth = (proj.abp)./(proj.Sp)
a = quantile(wealth,[0.1])
b = quantile(wealth,[0.2])
c = quantile(wealth,[0.5])
d = quantile(wealth,[0.9])

omega = weights_bin

indcc=findall(x-> x<= a[1], wealth ) #this is the index
P1 = mean(omega[indcc])
indcc2=findall(x-> a[1]<x<=b[1],wealth ) #this is the index
P2 = mean(omega[indcc2])
indcc3=findall(x-> b[1]<x<=c[1],wealth ) #this is the index
P3 = mean(omega[indcc3])
indcc4=findall(x-> c[1]<x<=d[1],wealth ) #this is the index
P4 = mean(omega[indcc4])
indcc5=findall(x-> d[1]<x ,wealth ) #this is the index
P5 = mean(omega[indcc5])

P = [P1, P2, P3, P4, P5]
t = 1:length(P)
x = skipmissing(P)
mean(skipmissing(P))

scatter(P, label="")
nonmissing = isfinite.(t .+ P)
plot(t[nonmissing], P[nonmissing], label="", linewidth = 2, dpi=300)
xlabel!("Wealth level")
ylabel!("Average Pareto Weights")
ylims!(0,4)
#savefig("weights_france_wealth_average.png")


## Alternative way to compute the Pareto weights 
M8= zeros(proj.Nbin,ns)
for i=1:proj.Nbin
    M8[i,proj.ytype[i]]=1
end

weights_bin_y= M8*weights

plot(wealth./1000,states[proj.ytype],weights_bin_y,st=:surface, legend =(0.1,0.1), zlabel=("Pareto weights"), zguidefontrotation=90, yguidefontrotation=-18, xguidefontrotation=50, dpi=300 ,camera=(-60,30))
zlims!(0.0,5.0)
xlabel!("Wealth in k of €")
ylabel!("Productivity level")

plot(wealth./1000,wage./1000,weights_bin_y,st=:surface, legend =(0.1,0.1), zlabel=("Pareto weights"), zguidefontrotation=90, yguidefontrotation=-18, xguidefontrotation=50, dpi=300 , camera=(-60,30))
zlims!(0.0,5.0)
xlabel!("Wealth in k of €")
ylabel!("Average wage per year in k of €")

omega_y = weights_bin_y

indcc=findall(x-> x<= a[1], wealth ) #this is the index
P1 = mean(omega_y[indcc])
indcc2=findall(x-> a[1]<x<=b[1],wealth ) #this is the index
P2 = mean(omega_y[indcc2])
indcc3=findall(x-> b[1]<x<=c[1],wealth ) #this is the index
P3 = mean(omega_y[indcc3])
indcc4=findall(x-> c[1]<x<=d[1],wealth ) #this is the index
P4 = mean(omega_y[indcc4])
indcc5=findall(x-> d[1]<x ,wealth ) #this is the index
P5 = mean(omega_y[indcc5])

P = [P1, P2, P3, P4, P5]
t = 1:length(P)
x = skipmissing(P)
mean(skipmissing(P))

scatter(P, label="")
nonmissing = isfinite.(t .+ P)
plot(t[nonmissing], P[nonmissing], label="", linewidth = 2, dpi=300)
xlabel!("Wealth level")
ylabel!("Average Pareto Weights")
ylims!(0,4)
#savefig("weights_france_wealth_average.png")


## Plotting the Weights by History
plot(weights_bin, legend=:topleft, line = (1,:dash, [:blue]), label = "Weights for each history", dpi=300)
#plot!(weights_bin_y,legend=:topleft, line = (1,:solid, [:red]), label = "Weights for productivity", dpi=300)
ylims!(-0.2,5)
xlabel!("Histories")
ylabel!("Pareto weights")
savefig("weights_france_comparison.png")


weights_average_bin =weights_bin


## Initialize an array to hold the sums for each column
sums = Float64[]

## Iterate over each column
for col in 1:size(M8, 2)
    current_sum = 0.0
    for row in 1:size(M8, 1)
        if M8[row, col] == 1
            current_sum += weights_average_bin[row]
        elseif M8[row, col] == 0 && current_sum > 0
            # Store the current sum when hitting a zero after summing ones
            push!(sums, current_sum)
            current_sum = 0.0  # Reset the sum for the next group
        end
    end
    # If we finished with a non-zero sum, add it
    if current_sum > 0
        push!(sums, current_sum)
    end
end

## Output the result
println(sums)

weighted_values = weights_average_bin .* proj.Sp


## Initialize an array to hold the normalized sums for each column
normalized_sums = Float64[]


## Iterate over each column in M8
for col in 1:size(M8, 2)
    current_sum = 0.0
    current_proj_sum = 0.0  #sum of proj.Sp for normalization
    for row in 1:size(M8, 1)
        if M8[row, col] == 1
            #add to current_sum and current_proj_sum if M8[row, col] is 1
            current_sum += weighted_values[row]
            current_proj_sum += proj.Sp[row]
        elseif M8[row, col] == 0 && current_sum > 0
            #store the normalized sum when hitting a zero
            push!(normalized_sums, current_sum / current_proj_sum)
            current_sum = 0.0  #reset the sum for the next group
            current_proj_sum = 0.0  #reset proj.Sp sum for normalization
        end
    end
    #if there’s a remaining sum, add the normalized sum
    if current_sum > 0
        push!(normalized_sums, current_sum / current_proj_sum)
    end
end


## Output the result
println(normalized_sums)

normalized_sums = normalized_sums./mean(normalized_sums)


## Plotting the Comparison between the Weights by History and the Non-Parametri Weightss
plot(weights,line = (2,:solid, [:blue]), label = "Weights")
plot!(normalized_sums, line = (2, :dash, :red), label = "Average Weights by History")
plot!(10*disty,legend=:topleft, line = (2, :dash, [:darkgreen]),fill = (0, :lightgreen),fillalpha = 0.3,dpi = 300, label = "Distribution")
ylims!(0,5)
savefig("weights_france_comparison_history.png")
