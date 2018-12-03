# Author:  Marisol Monterrubio-Velasco
# Contact: marisol.monterrubio@bsc.es
# GNU GENERAL PUBLIC LICENSE
# ==========================
# Version 3, 29 June 2007
# Copyright © 2007 Free Software Foundation, Inc. <http://fsf.org/>

function contaravalan_Asp(nx,ny,VectorP,VecAsperi,VecPhi,rho,k)

# """
# contaravalan_Asp.jl : Compute the number of cells that overpass the threshold load value (STAGE 1). But also in this version the cell that fails is choosen considering not only its load but also the strength criterion defined in TREMOL algorithm (STAGE 2)

#Parameters
# - `nx::Integer`: number of cells in X-axis of the domain  \Omega
# - `ny::Integer`: number of cells in Y-axis of the domain \Omega 
# - `VectorP::Array`: size(nbox_x,nbox_y), Matrix of the load values. Dynamical matrix because changes the value of some cells at each step. 
# - `VecAsperi::Array`: size(nbox_x,nbox_y), Matrix of the strength value.  This matrix evolves at each step 
# - `VecPhi:Array`: size(nbox_x,nbox_y), Matrix of the load-tranfer values.
# - `rho::Integer`: Weibull exponent, we take it as a constant 30
# - `k::Integer`: step  of number

# return
# - `contador::Integer`: number of cells that overpass the threshold load value
# - `Nflag1::Integer` : code that indicates if the rupture is normal or avalanche 
# - `maxtempo1::Float` : maximum load value
# - `iout1::Integer`: x coordinate in the array of the chosen cell to fail
# - `jout1::Integer`: y coordinate in the array of the chosen cell to fail
# - `VecAsperi::Array`: updated matrix of the strength
# - `iout::Integer`: x coordinate in the array of the chosen cell to fail at STAGE 1
# - `jout::Integer`: Y coordinate in the array of the chosen cell to fail at STAGE 1
# """

    contador=0
    Nflag=0
    paramVecAval=0
    VecAval=zeros(nx*ny)
    VecFueraAsperi = zeros(nx*ny)
    iout = 0
    jout = 0
    maxtempo = 0.0
    maxtempo1 = 0.0
    iout1 = 0
    jout1 = 0

# STAGE 1: in this loop the number of cells that overpass the threshold load value are counted.  
    
    VecPNotBord = VectorP[2:nx,2:ny]
    VectorPArray = VecPNotBord[:]  
    
    VecAvalIndex = find(map(x -> x >= 1.0, VectorPArray)) 
    VecAval = VectorPArray[VecAvalIndex[:]]
    contador = length(VecAvalIndex)
      
    if contador == 0
        Nflag = 0
    else
        Nflag = 1
    end

    Nflag1 = 0

    if contador >= 1  # If one or more cells overpass the threshold load we go to this if statment (Avalanche event)o avalancha
        paramVecAval=maximum(VecAval[1:contador])
        B = ind2sub(size(VectorP),indmax(VectorP))
        iout = getindex(B,1)
        jout = getindex(B,2)      
        maxtempo = VectorP[iout,jout]
        iout1 = 0
        jout1 = 0
        VecAsperiNotBord = VecAsperi[2:nx,2:ny]
        VectorAsperiArray = VecAsperiNotBord[:]
        
# The strength of the cell is checked      

# STAGE 2: If the cell has a strength larger of 1, the stage two starts and we search other cell that has a strength equal to one and a larger load to fail.
       
        if VecAsperi[iout,jout] > 1  # If the chosen element in STAGE 1 has a strength larger than one, then that cell is weakened 
            VecAsperi[iout,jout] = VecAsperi[iout,jout] - 1         
            contFuera = 0      
            VectorPcopy = copy(VectorP)
            VectorPcopy[iout,jout] = -12 
            VecPNotBordcopy = VectorPcopy[2:nx,2:ny]
            VectorPArraycopy = VecPNotBordcopy[:]  
      
            for i = 1:length(VectorAsperiArray)           
                if  VectorAsperiArray[i] == 1 && VectorPArraycopy[i] > 0.0  # search among all the elements of the matrix that have stregth = 1 not include the element chosen in STAGE 1
                    contFuera += 1
                    VecFueraAsperi[contFuera] = VectorPArraycopy[i]
                end
            end             
            if contFuera > 0        
                maxtempo1 =maximum(VecFueraAsperi[1:contFuera])
                if maxtempo1 >= 1.0  # if the item being searched has a load greater than 1
                    Nflag1=1    # flag to indicate if we choose an avalanche event    
                    for j = 2:ny
                        for i = 2:nx
                            if maxtempo1 == VectorP[i,j]
                                iout1 = i
                                jout1 = j
                                break
                            end
                        end 
                    end         
                end       
            else
                contador = -2
                iout1 = 0
                jout1 = 0
                Nflag1 = 0 
                maxtempo1 = 0
            end         
            if Nflag1 == 0 # a normal event with strength equal to one is choosen  
                r=0
                l=0
                sumarho=sum(VecFueraAsperi[1:contFuera].^rho)
                tempsmig2 = 1.0/sumarho
                probkFuera=zeros(contFuera)
                limiteinferior=0
                probab=rand()   # random number (0,1) 
                VecFueraAsperiRand = shuffle(VecFueraAsperi[1:contFuera])                                                  
                contador2 = 0
                eleminfalla = 0
                for i=1:contFuera
                    limiteinferior = contador2
                    probkFuera[i] = (VecFueraAsperiRand[i].^rho)*tempsmig2  # the probability pk is computed for each element 
                    contador2 += probkFuera[i]	  # a cumulative probability is increasing
                    if limiteinferior < probab < contador2    # a comparison between the random number `probab` and the cumulative probability is carry out
                        eleminfalla = VecFueraAsperiRand[i]
                        numOrdEleMin = i
                        break
                    end
                end    
                for j=2:ny
                    for i=2:nx
                        if  eleminfalla == VectorP[i,j]  # the chosen element is find in the load array
                            fhi=VecPhi[i,j]
                            iout1 = i
                            jout1 = j
                            maxtempo1 = VectorP[i,j]
                            break
                        end
                    end
                end
            end              
        elseif VecAsperi[iout,jout] == 1          
            iout1 = iout
            jout1 = jout
            maxtempo1 = maxtempo
            Nflag1 = 1
        end
    elseif contador == 0  # normal event   
        rN=0
        lN=0
        suma = 0
        sumarho = 0
        tempsmig = 0
        vecNorm = zeros(length(VectorPArray))
        contNorm = 0        
        for i = 1:length(VectorPArray)
            if VectorPArray[i] != -1.0 
                contNorm += 1
                suma += VectorPArray[i]
                sumarho += (VectorPArray[i]^rho)
                vecNorm[contNorm] = VectorPArray[i]
            end    
        end  
        if contNorm == 0
            contador = -2           
            iout1 = 0
            jout1 = 0
            Nflag1 = 0 
            maxtempo1 = 0          
        elseif contNorm >= 0    
            tempsmig = 1.0/sumarho
            probkFuera = zeros(nx*ny)
            limiteinferior = 0
            probab = rand()      
            contador2 = 0
            ioutN = 0
            joutN = 0
            contFueraX = 0
            eleminfalla = 0      
            vecNormRand = shuffle(vecNorm[1:contNorm])    
            for i = 1:contNorm
                limiteinferior = contador2
                probkFuera[i] = (vecNormRand[i].^rho)*tempsmig   # the probability pk is computed for each element 
                contador2 += probkFuera[i]	       # a cumulative probability is increasing
                if limiteinferior < probab < contador2   # a comparison between the random number `probab` and the cumulative probability is carry out
                    eleminfalla = vecNormRand[i]
                    break
                end        
            end   
            for j = 2:ny
                for i = 2:nx
                    if  VectorP[i,j] == eleminfalla   # the chosen element is find in the load array
                        fhi=VecPhi[i,j]
                        ioutN = i
                        joutN = j
                        maxtempo1 = VectorP[i,j]
                    end
                end
            end
        
            contFuera = 0
            iout = ioutN 
            jout = joutN  
            VecFueraAsperi2 = zeros(nx*ny)     

#   if the chosen event has a strength larger than one VecAsperi>1          
            if VecAsperi[iout,jout] > 1
                maxtempo0 = 0
                VecAsperi[iout,jout] = VecAsperi[iout,jout] - 1        
                for j = 2:ny
                    for i = 2:nx
                        if VecAsperi[i,j] == 1 && i != iout && j != jout
                            if VectorP[i,j] > 0.0
                                contFuera += 1
                                VecFueraAsperi2[contFuera] = VectorP[i,j]
                            end
                        end
                    end
                end             
                probkFuera=zeros(contFuera)
                sumarho2=sum(VecFueraAsperi2[1:contFuera].^rho)
                tempsmig2 = 1.0/sumarho2
                contador3 = 0
                eleminfallaN=0
                limiteinferior=0
                iout0 = 0
                jout0 = 0
                probab = rand()        
                VecFueraAsperi2Rand = shuffle(VecFueraAsperi2[1:contFuera])       
                for i=1:contFuera
                    limiteinferior = contador3              
                    probkFuera[i] = (VecFueraAsperi2Rand[i].^rho)*tempsmig2   
                    contador3 += probkFuera[i]	        
                    if limiteinferior < probab < contador3    
                        eleminfallaN = VecFueraAsperi2Rand[i]
                        break
                    end
                end
                for j=2:ny
                    for i=2:nx
                        if  eleminfallaN == VectorP[i,j]    
                            fhi=VecPhi[i,j]
                            iout1= i
                            jout1 = j
                            maxtempo1 = VectorP[i,j]
                            VecAsperi[iout1,jout1] = 1
                            break
                        end
                    end
                end
                Nflag1 = 0  
            elseif  VecAsperi[iout,jout] == 1 
                iout1 = iout
                jout1 = jout
                maxtempo1 = maxtempo1
                Nflag1 = 0 
            end
        end
    end

    return contador,Nflag1,maxtempo1,iout1,jout1,VecAsperi,iout,jout
end

