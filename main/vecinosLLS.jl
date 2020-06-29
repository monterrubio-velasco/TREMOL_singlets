# Author:  Marisol Monterrubio-Velasco
# Contact: marisol.monterrubio@bsc.es
# GNU GENERAL PUBLIC LICENSE
# ==========================
# Version 3, 29 June 2007
# Copyright © 2007 Free Software Foundation, Inc. <http://fsf.org/>

include("veciProhiDiag.jl")

function vecinosLLS(vector,vecEsfTot)

# """
# vecinosLLS.jl : This function distributed the load for the orthogonal neigbors (N,S,E,W) that are allowed to received load. 
# 
#Import function
# - include(joinpath("/YOURPATH/TREMOL_singlets/main/","veciProhiDiag.jl"))

#Parameters
# - `vector::Array`:vector of nine positions that contains the neighbors and failed cell, being:   vector[9], vector[3], vector[7] and vector[1] diagonal neighbors; vector[8], vector[2], vector[4] and vector[6] perpendicular neighbors.
# - `vecEsfTot::Float`: Load value of the failed cell.    

# return
# - `vectclon::Array`: Updated vector of nine positions after the load transfer of the failed cell.
# """

    vecEsfDiag = vecEsfTot*0.02  # 2% of load is distributed to the neighbors located at NE-NW-SE-SW
    factDiag = vecEsfDiag/4      # divided the distributed load by 4 neighbors NE-NW-SE-SW
    vecEsf = vecEsfTot*0.98      # 98% of load is distributed to the neighbors located at N-S-E-W
    factPerpe = vecEsf/4.0       # divided the distributed load by 4 neighbors N-S-E-W

# If statment to distributed the load only to the no prohibited neighbors 

#  Orthogonal neighbors

# Option 1: all the N-S-E-W are forbidden
    if vector[2] < 0 && vector[4]< 0 && vector[6] < 0 && vector[8] < 0
        vector[5] = -1.0
    end
# Option 2: Three neighbors are forbidden
    if vector[4] < 0 && vector[6] < 0 && vector[8] < 0 && vector[2] > 0.0
        vector[2] = vector[2] + vecEsf    
    elseif vector[6] < 0 && vector[8] < 0 && vector[2] < 0 && vector[4] > 0.0
        vector[4] = vector[4] + vecEsf      
    elseif vector[8] < 0 && vector[2] < 0 && vector[4] < 0 && vector[6] > 0.0
        vector[6] = vector[6] + vecEsf        
    elseif vector[2] < 0 && vector[4] < 0 && vector[6] < 0 && vector[8] > 0.0
        vector[8] = vector[8] + vecEsf
    end
# Option 3: Two neighbors are forbidden
    if vector[8] < 0 && vector[2]< 0 && vector[4] > 0.0 && vector[6] > 0.0
        vector[4] = vector[4] + (vecEsf/2)
        vector[6] = vector[6] + (vecEsf/2)    
    elseif vector[2]< 0 && vector[4] < 0 && vector[8] > 0.0 && vector[6] > 0.0
        vector[8] = vector[8] +(vecEsf/2)
        vector[6] = vector[6] + (vecEsf/2)      
    elseif vector[2]< 0 && vector[6] < 0 && vector[8] > 0.0 && vector[4] > 0.0
        vector[8] = vector[8] + (vecEsf/2)
        vector[4] = vector[4] + (vecEsf/2)      
    elseif vector[4] < 0 && vector[6] < 0 && vector[8] > 0.0 && vector[2]> 0.0
        vector[8] = vector[8] + (vecEsf/2)
        vector[2]= vector[2]+ (vecEsf/2)      
    elseif vector[8] < 0 && vector[4] < 0 && vector[2]> 0.0 && vector[6] > 0.0
        vector[2]= vector[2]+ (vecEsf/2)
        vector[6] = vector[6] + (vecEsf/2)      
    elseif vector[8] < 0 && vector[6] < 0 && vector[2] > 0.0 && vector[4] > 0.0
        vector[2] = vector[2]+(vecEsf/2)
        vector[4] = vector[4]+(vecEsf/2)      
    end
# Option 4: only one neighbor is forbidden 
    if vector[8] < 0 && vector[2] > 0.0 && vector[4] > 0.0 && vector[6] > 0.0
        vector[2]= vector[2] + (vecEsf/3)
        vector[4] = vector[4] +  (vecEsf/3)
        vector[6] = vector[6] +  (vecEsf/3)      
    elseif vector[2] < 0 && vector[8] > 0.0 && vector[4] > 0.0 && vector[6] > 0.0
        vector[8] = vector[8] +  (vecEsf/3)
        vector[4] = vector[4] + (vecEsf/3)
        vector[6] = vector[6] + (vecEsf/3)      
    elseif vector[6] < 0 && vector[2] > 0.0 && vector[8] > 0.0 && vector[4] > 0.0
        vector[2]= vector[2]+  (vecEsf/3)
        vector[8] = vector[8] + (vecEsf/3)
        vector[4] = vector[4] +  (vecEsf/3)      
    elseif vector[4] < 0 && vector[2] > 0.0 && vector[8] > 0.0 && vector[6] > 0.0
        vector[2] = vector[2]+  (vecEsf/3)
        vector[8] = vector[8] +  (vecEsf/3)
        vector[6] = vector[6] +  (vecEsf/3)
    end

#  Diagonal neighbors
    D=zeros(4)
    D[1]=vector[9]
    D[2]=vector[7]
    D[3]=vector[1]
    D[4]=vector[3]
    B = veciProhiDiag(vector,vecEsfDiag,D)
    vector[9]=B[1]
    vector[7]=B[2]
    vector[1]=B[3]
    vector[3]=B[4]
    vectclon=copy(vector)
    
    return vectclon
end
