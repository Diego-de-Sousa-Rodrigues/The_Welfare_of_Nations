function Weight_Matrix(Modeling::AiyagariModelEGM,WeightP::Array{T},WeightS::Array{T})where{T <: Float64}

    ny = Modeling.ns       # number of productivity levels
    piy = Modeling.disty   # shares by levels

  
 # normalization with political weights, for the weighted sum to be 1
    #WeightS = WeightS/sum(WeightS.*piy.*WeightP)

    WM = zeros(Float64,ny,ny) # matrix with political weights as a loading factor
    #Diagonal(1.0 ./(WeightP.*piy))
    #diagv = 1.0 ./(WeightP.*piy)
    WM  = Diagonal(vec(1.0 ./(WeightP.*piy))) + (WeightS - ones(length(WeightS),1) )*ones(1,length(WeightS))

   
    C2 = sum(piy.*WeightP.*WeightP)
    
    WM2 = zeros(Float64,ny,ny) # Good one: matrix without political weights as a loading factor
    for i=1:ny
        for j=1:ny
            WM2[i,j] = WeightP[i]*piy[i]*((i==j)*1/(WeightP[i]*piy[i]) + (WeightS[j]-1)*WeightP[i]/C2 )
        end
    end

    WM2s = zeros(Float64,ny,ny) # different weights
    for i=1:ny
        for j=1:ny
            WM2s[i,j] = (WeightS[j]-1)*WeightP[i]/C2 
        end
    end

    return WM,WM2,WM2s
    end