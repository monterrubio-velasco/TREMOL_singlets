# Author:  Marisol Monterrubio-Velasco
# Contact: marisol.monterrubio@bsc.es

include("vecinosLLS.jl")


function distmayorLLS_Asp(vector,fhi)

# """
# distmayorLLS_Asp.jl : Function that distribute the load following the avalanche events algorithm (avalanche event is the cell that overpass their threshold values)

#Parameters
# - `vector::Array`:vector of nine positions that contains the neighbors and failed cell, being:   vector[9], vector[3], vector[7] and vector[1] diagonal neighbors; vector[8], vector[2], vector[4] and vector[6] perpendicular neighbors. Finally vector[5] is the failed cell.     
# - `fhi::Float`: load-transfer value assigned to the failed cell. 

#return
# - `vector::Array`: updated vector of nine positions after the load transfer of the failed cell.
# """

    vecEsfTot = vector[5]*fhi   #distributed load from failed cell to their neighbors
    vecEsfDiag = vecEsfTot*0.02 # 2% of load is distributed to the neighbors located at NE-NW-SE-SW
    factDiag = vecEsfDiag/4    # divided the distributed load by 4 neighbors NE-NW-SE-SW
    vecEsf = vecEsfTot*0.98    # 98% of load is distributed to the neighbors located at N-S-E-W
    factPerpe=vecEsf/4         # divided the distributed load by 4 neighbors N-S-E-W
    nflag1 = 0 
# Checking if some neighbor is forbidden
    if vector[1] < 0 || vector[2] < 0 || vector[3] < 0 || vector[4] < 0 || vector[6] < 0 || vector[7] < 0 || vector[8] < 0 || vector[9] < 0
        vectclon = vecinosLLS(vector,vecEsfTot)
        vectclon[5] = -1.0
        nflag1=1
    else	    

        vector[5] = -1.0
# neighbors orthogonal located NE-NW-SE-SW
        vector[8] += factPerpe
        vector[2] += factPerpe
        vector[4] += factPerpe
        vector[6] += factPerpe         
# neighbors diagonally located NE-NW-SE-SW
        vector[9] += factDiag
        vector[3] += factDiag
        vector[7] += factDiag
        vector[1] += factDiag       	              	      
    end
    if nflag1==1
        vector=copy(vectclon)
    end       
    return vector

end

