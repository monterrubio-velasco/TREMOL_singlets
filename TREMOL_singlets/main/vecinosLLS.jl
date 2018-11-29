# Author:  Marisol Monterrubio-Velasco
# Contact: marisol.monterrubio@bsc.es

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

    vecEsfDiag = vecEsfTot*0.02  # 2% de esfuerzo a las vecinas NE-NW-SE-SW
    factDiag = vecEsfDiag/4      #las cuatro direcciones

  # Se reparte en N-S-E-W un porcentaje mayor de Esfuerzo
    vecEsf = vecEsfTot*0.98  # 98% de esfuerzo a las vecinas N-S-E-W
    factPerpe = vecEsf/4.0

# Bucle que mira si los vecinos son prohibidos. Si lo son el esfuerzo
# se distribuye entre los vecinos no prohibidos, segun el número de vecinos a repartir

  # Opción 1: TODOS LOS VECINOS N-S-E-W SON PROHIBIDOS
#   if vector[2] < 0 && vector[4]< 0 && vector[6] < 0 && vector[8] < 0
#  	   vector[5] = -1.0
# 	end
 # Opcion 2: TRES VECINOS SON PROHIBIDOS
    if vector[4] < 0 && vector[6] < 0 && vector[8] < 0 && vector[2] > 0.0
        vector[2] = vector[2] + vecEsf    
    elseif vector[6] < 0 && vector[8] < 0 && vector[2] < 0 && vector[4] > 0.0
        vector[4] = vector[4] + vecEsf      
    elseif vector[8] < 0 && vector[2] < 0 && vector[4] < 0 && vector[6] > 0.0
        vector[6] = vector[6] + vecEsf        
    elseif vector[2] < 0 && vector[4] < 0 && vector[6] < 0 && vector[8] > 0.0
        vector[8] = vector[8] + vecEsf
    end
	# Opcion 3: DOS VECINOS SON PROHIBIDOS	
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
  # Opcion 4: UN VECINO ES PROHIBIDO
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

   # Los vecinos diagonales
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
