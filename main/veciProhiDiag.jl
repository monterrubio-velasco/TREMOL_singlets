# Author:  Marisol Monterrubio-Velasco
# Contact: marisol.monterrubio@bsc.es
# GNU GENERAL PUBLIC LICENSE
# ==========================
# Version 3, 29 June 2007
# Copyright © 2007 Free Software Foundation, Inc. <http://fsf.org/>

function veciProhiDiag(vector,vecEsfDiag,A)

# """
# veciProhiDiag.jl : Compute the number of cells that overpass the threshold load value. But also in this version it choose the cells that is ruptured because the strength criterion defined in TREMOL algorithm

#Parameters
# - `vector::Array`:vector of nine positions that contains the neighbors and failed cell, being:   vector[9], vector[3], vector[7] and vector[1] diagonal neighbors; vector[8], vector[2], vector[4] and vector[6] perpendicular neighbors.`
# - `vecEsfDiag::Float`: Load value will be transfer to the Diagonal neighbors. `
# - `A::Array`:vector of four positions that contains the diagonal neighbors

#return
# - `B::Array`: updated vector of four positions that contains the new amount of load given to the diagonal neighbors
# """

    contador=0
    for i=1:length(A)
        if A[i] < 0
            contador=contador+1
        end
    end
    if contador != length(A)
        num=length(A)-contador
        newsigma=vecEsfDiag/num
        for  i=1:length(A)
            if A[i] > 0
                A[i] = A[i] + newsigma
            end
        end    
    end
    B=A
    return B
end
