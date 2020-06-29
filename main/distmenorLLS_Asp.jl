# Author:  Marisol Monterrubio-Velasco
# Contact: marisol.monterrubio@bsc.es
# GNU GENERAL PUBLIC LICENSE
# ==========================
# Version 3, 29 June 2007
# Copyright © 2007 Free Software Foundation, Inc. <http://fsf.org/>

include("vecinosLLS.jl")

function distmenorLLS_Asp(vector1,fhi)  

# """
# distmenorLLS_Asp.jl : Function that distributes the load following the normal-events algorithm (when any cell ooverpass their threshold values)

#Parameters
# - `vector1::Array`:vector of nine positions that contains the neighbors and failed cell, being:   vector1[9], vector1[3], vector1[7] and vector1[1] diagonal neighbors; vector1[8], vector1[2], vector1[4] and vector1[6] perpendicular neighbors. Finally vector1[5] is the failed cell.     
# - `fhi::Float`: load-transfer value assigned to the failed cell 

# return
# - `vecreload::Array`: updated vector of nine positions after the load transfer of the failed cell.
# """
    vecEsfTot=vector1[5]*fhi
    vecEsfDiag=vecEsfTot*0.02   # 2% of load is distributed to the neighbors located at NE-NW-SE-SW
    factDiag=vecEsfDiag/4    # divided the distributed load by 4 neighbors NE-NW-SE-SW
    vecEsf=vecEsfTot*0.98   # 98% of load is distributed to the neighbors located at N-S-E-W
    factPerpe=vecEsf/4      # divided the distributed load by 4 neighbors N-S-E-W
    nflag1=0  
    vector1[5] = 0.0

# Checking if some neighbor is forbidden
    if vector1[1] < 0 || vector1[2] < 0 || vector1[3] < 0 || vector1[4] < 0 || vector1[6] < 0 || vector1[7] < 0|| vector1[8] < 0|| vector1[9] < 0
        vectclon = vecinosLLS(vector1,vecEsfTot)
        nflag1=1
    else
      vector1[8] += factPerpe
      vector1[2] += factPerpe
      vector1[4] += factPerpe
      vector1[6] += factPerpe
# neighbors diagonally located
      vector1[9] += factDiag
      vector1[3] += factDiag
      vector1[7] += factDiag
      vector1[1] += factDiag
    end
    if nflag1==1
      vector1=copy(vectclon)
      vector1[5] = 0.0
    end
    vecreload=copy(vector1)
    
    return vecreload
end
