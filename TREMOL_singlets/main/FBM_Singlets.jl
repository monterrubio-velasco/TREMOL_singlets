# Author:  Marisol Monterrubio-Velasco
# Contact: marisol.monterrubio@bsc.es
# GNU GENERAL PUBLIC LICENSE
# ==========================
# Version 3, 29 June 2007
# Copyright © 2007 Free Software Foundation, Inc. <http://fsf.org/>

include("contaravalan_Asp.jl")
include("distmayorLLS_Asp.jl") 
include("distmenorLLS_Asp.jl")   

function FBM_Singlets(nbox_x,nbox_y,smin,VecPosi,MatrizStrengtInitial,VecAsperi,fhiFuera,fhiDentro,VecPhi)

# """
# FBM_Singlets.jl : This function carry out the FBM asperity algorithm 
# 
#Import function
# - include(joinpath("/YOURPATH/TREMOL_singlets/main/","contaravalan_Asp.jl"))
# - include(joinpath("/YOURPATH/TREMOL_singlets/main/","distmayorLLS_Asp.jl"))  
# - include(joinpath("/YOURPATH/TREMOL_singlets/main/","distmenorLLS_Asp.jl"))
# 
#Parameters
# - `nbox_x::Integer`: number of cells in X-axis of the domain  \Omega
# - `nbox_y::Integer`: number of cells in Y-axis of the domain \Omega 
# - `smin::Integer`: number of steps that realize the algorithm 
# - `VecPosi::Array`: size(nbox_x,nbox_y), Matrix of the load values. Dynamical matrix because changes the value of some cells at each step. 
# - `MatrizStrengtInitial::Array`: size(nbox_x,nbox_y), Initial Matrix of the strength value.
# - `VecAsperi::Array`: size(nbox_x,nbox_y), Matrix of the strength value.  This matrix evolves at each step 
# - `fhiFuera::Float`: load-transfer value assigned to the cells located in the domain
# - `fhiDentro::Float`: load-transfer value assigned to the cells located in the asperity domain  
# - `VecPhi::Array`: size(nbox_x,nbox_y), Matrix of the load-tranfer values.
# 
# return
# - `vectk1::Array`: size(smin+10,12). Raw data that contains the results of the transfer, accumulation and rupture process following the FBM rules. At each step   each row contains: 
#     1. k:number of step
#     2. acumt: cumulative time (T_k = sum(tiempo[1:k]))
#     3. tiempo[k]: inter-event time [dimensionless]
#     4. (0 or 1): counter that indicates if is normal or avalanche event
#     5. suma: sum of the load in all the cells
#     6. sumarho= 1/tiempo[k]  (rupture rate)
#     7. parametrosigma: load value of the cell chosen to fail
#     8. a: coordinate in the X-axis, of the cell chosen to fail
#     9. b: coordinate in the Y-axis, of the cell chosen to fail
#     10. fhi: load-transfer value
#     11-12. (iorigi,jorigi): coordinates of the cell chosen to fail at the first stage of the searching algortithm defined in the function contaravalan_Asp. 
#     If their stregth is equal to one it fails if is larger the searching algorithm will chose a cell taken into account the cells with strength equal to one. 
# - `vectparamestad::Array`:  size(smin+10,2). Mean and Standard deviation of the load in the system, computed at each time step considering only the cells active 
# - `VecPosiFinal::Array`: size(nbox_x,nbox_y). Array of the final configuration of VecPosi.
# """

#   the border is added to the initial arrays VecAsperi,VecPhi, VecPosi 
    srand(1)
    Z=zeros(nbox_x, 1)
    A=[Z VecPosi Z]
    Z1=zeros(1, length(A[1,:]))
    C=[Z1; A; Z1]  
    
    A1=[Z VecAsperi Z]
    Z2=zeros(1, length(A1[1,:]))
    C1=[Z2; A1; Z2]
  
    A2=[Z VecPhi Z]
    Z3=zeros(1, length(A2[1,:]))
    C2=[Z3; A2; Z3]
  
    A3 = [Z MatrizStrengtInitial Z]
    Z4 = zeros(1, length(A3[1,:]))
    C3 = [Z4; A3; Z4]
  
    VecPosi = copy(C)
    VecAsperi = copy(C1)
    VecPhi = copy(C2)
    VecMatrizStrenIniti = copy(C3)

    Nboxy=nbox_y+1
    Nboxx=nbox_x+1

#   All vectors are initialized at zero
    tiempo = zeros(smin+10)
    VecFlag = zeros(smin+10)
    vectk1 = zeros(smin+10,12)
    vectFinSimul = zeros(smin,2)
    vectRepSimul = zeros(smin,3)
    vectparamestad = zeros(smin+1,2)
    vectparamestadFuera = zeros(smin+1,2)
    vectparamestadDentro = zeros(smin+1,2)
    VecDentroLoad = zeros(Nboxy*Nboxx)
    VecFueraLoad =  zeros(Nboxy*Nboxx) 
    contFin = 0
    contRep = 0
    label = 0
    contini = 0
    nflag = 0   # counter that Contador que indica si es EventoNormal (NFLAG=0) o EventoAvalancha (NFLAG=1)
    flagindica = 0
    contava = 0
    contAval = 0
    contNumAval = 0
    conGRAva = 0
    numCodntMag = 0
    numAvalMag = 0
    
#   Weibull coefficient 
    rho=30      
    TempMean = 0
    VecMeanTempo = zeros(Nboxx*Nboxy)
    contadorDentro = 0
    contadorFuera = 0  
    tamatresf = (Nboxx-2)*(Nboxy-2)
    sumaexponente=sum(VecPosi[2:Nboxx,2:Nboxy])
    sumapotencia=sum(VecPosi[2:Nboxx,2:Nboxy].^rho) 
    acumt=0.0 # Cumulative time. The initial time is t=0, acumt=0.
    numcontador=0  
        
#   k step counter initialized at zero  
    k = 1  
    vectparamestad[k,1] = mean(VecPosi[2:Nboxx,2:Nboxy])
    vectparamestad[k,2] = std(VecPosi[2:Nboxx,2:Nboxy])
    tiempo[k] = 1.0/sumapotencia  # interevent time
    VecFlag[k] = 0
    acumt += tiempo[k] # cumulative time 
  
  ####||||||||||||||||||||||||||||######################
      

#   first row of the synthetic raw database 
    vectk1[k,1] = k 
    vectk1[k,2] = acumt
    vectk1[k,3] = tiempo[k]
    vectk1[k,4] = VecFlag[k]
    vectk1[k,5] = sumaexponente
    vectk1[k,6] = sumapotencia   
    
    indicenorm=0
    contNorm=1
# Vector que guarda los eventos Normales
  
    nflagUltimo = 0   
    a=0
    b=0
    contNucIni = 1

########## Starts the algorithmic loop ################
    while k < smin  
        k = k+1  # step increases
        println(k)
        
        contador,nflag,maxLoad,a,b,VecAsperiOut,iorigi,jorigi = contaravalan_Asp(Nboxx,Nboxy,VecPosi,VecAsperi,VecPhi,rho,k)    # function contaravalan_Asp
        VecAsperi = VecAsperiOut 
      
        if contador == -2  # sistema agotado
            smin = k-1
            nflagUltimo = nflag
        end
      
        if a == 0 || b == 0  # sistema agotado
            k = k-1
            println("a = 0 or b = 0")
            contador = -2
            pause()
        end
      
        if contador >= 1  #Se está en un evento Avalancha

            if nflag == 1

                parametrosigma = maxLoad # La celda que supero su esfuerzo umbral es la elegida para fallar
                FailCell = VecPosi[a,b]
                contaBucleAval = 0
                NeighbordCell = zeros(9)
            
                for jk = b - 1 : b + 1 
                    for ik = a - 1 : a + 1
                        contaBucleAval += 1
                        NeighbordCell[contaBucleAval] = VecPosi[ik,jk]                                                                
                    end	
                end  
        
                fhi = VecPhi[a,b]
                VecFlag[k] = nflag         
                vecnew = distmayorLLS_Asp(NeighbordCell,fhi)  
                            
                contaBucleAval = 0
                for jk = b - 1 : b + 1
                    for ik = a - 1 : a + 1 
                        contaBucleAval += 1
                        VecPosi[ik,jk] = vecnew[contaBucleAval]
                    end	
                end

                VecPosiArray = VecPosi[2:Nboxx,2:Nboxy]
                B = filter(x -> x != -1.0, VecPosiArray)
                suma = sum(B)
                sumarho = sum(B.^rho)   
                tiempo[k] = (1.0/sumarho)
                acumt += tiempo[k]
                vectk1[k,:] = [k,acumt,tiempo[k],1,suma,sumarho,parametrosigma,a,b,fhi,iorigi,jorigi]                    
                      
                vectparamestad[k,1] = mean(B)
                vectparamestad[k,2] = std(B)
                  
            elseif nflag == 0     # se elige un evento Normal con dureza = 1
      
                VecFlag[k] = nflag       
          
                if VecFlag[k-1] == 1  
                    for j = 2:Nboxy
                        for i = 2:Nboxx          
                            if VecMatrizStrenIniti[i,j] != 1 &&  VecPosi[i,j] == -1.0 
                                VecPosi[i,j] = -1.0   ##  cells that were avalanche events within the asperity remain not allow to receive any load                           
                            elseif VecMatrizStrenIniti[i,j] == 1 &&  VecPosi[i,j] == -1.0 
                                VecPosi[i,j] = 0.0   ##  cells that were avalanche events outside the asperity are reactivated with a load value equal to 0.0
                            end                          
                        end
                    end
                end
          
                NeighbordCell = zeros(9)
                contaBucleNorm = 0    
                for jk = b-1 : b+1
                    for ik = a-1 : a+1 
                        contaBucleNorm += 1
                        NeighbordCell[contaBucleNorm] = VecPosi[ik,jk]                                                                
                    end	
                end               

                fhi = VecPhi[a,b]  
                eleminfalla = NeighbordCell[5]           
                vecreload = distmenorLLS_Asp(NeighbordCell,fhi)        
                  
                contaBucleNorm = 0 
                for jk =  b - 1 : b + 1
                    for ik = a - 1 : a + 1
                        contaBucleNorm += 1
                        VecPosi[ik,jk] = vecreload[contaBucleNorm]           
                    end	
                end    
          
                suma = 0
                sumarho = 0
          
                B = filter(x -> x >= 0, VecPosi[2:Nboxx,2:Nboxy])
                suma = sum(B)
                sumarho = sum(B.^rho)
          
                tiempo[k] = 1.0/sumarho  
                acumt=acumt+tiempo[k]
                vectk1[k,:] = [k,acumt,tiempo[k],0,suma,sumarho,eleminfalla,a,b,fhi,iorigi,jorigi]
          
                contDisit = 0   
                
                if VecPhi[a,b] != fhiFuera
                    for ik = 1:length(vectFinSimul[:,1])
                        if a == vectFinSimul[ik,1] &&  b == vectFinSimul[ik,2]
                            break
                        else
                            contDisit +=1
                        end
                    end
                    if contDisit == length(vectFinSimul[:,1])
                        contFin +=1
                        if contFin > length(vectFinSimul[:,1])
                            break
                        else
                            vectFinSimul[contFin,1]=a
                            vectFinSimul[contFin,2]=b
                        end
                    end
                end
          
                TempMean = 0   
                VecMeanTempo = zeros(Nboxx*Nboxy)
                VecPosiArray = copy(VecPosi[2:Nboxx,2:Nboxy])
                VecPosiArray =  VecPosiArray[:]
          
                for i=1:length(VecPosiArray)
                    if VecPosiArray[i] >= 0.0
                        TempMean +=1
                        VecMeanTempo[TempMean]= VecPosiArray[i]
                    end
                end

              vectparamestad[k,1] = mean(VecMeanTempo[1:TempMean])
              vectparamestad[k,2] = std(VecMeanTempo[1:TempMean])                
              VecMeanTempo[1:TempMean]=0.0
              TempMean=0
              acumtref = acumt
              tiemporef = tiempo[k]
              indicenorm = indicenorm+1
              contNorm = contNorm+1
          end

        elseif contador == 0  
      
            VecFlag[k] = nflag                    
            contava=0    
    
            if VecFlag[k-1] == 1      
                for j = 2:Nboxy
                    for i = 2:Nboxx          
                        if VecMatrizStrenIniti[i,j] != 1 && VecPosi[i,j] == -1.0
                            VecPosi[i,j] = -1.0   ## cells that were avalanche events within the asperity remain not allow to receive any load 
                        elseif VecMatrizStrenIniti[i,j] == 1 && VecPosi[i,j] == -1.0
                            VecPosi[i,j] = 0.0   ## cells that were avalanche events outside the asperity are reactivated with a load value equal to 0.0
                        end         
                    end
                end
            end
        
            NeighbordCell = zeros(9)
            contaBucleNorm = 0 
            for jk = b - 1 : b + 1
                for ik =  a - 1 : a + 1
                    contaBucleNorm += 1
                    NeighbordCell[contaBucleNorm] = VecPosi[ik,jk]                                                                
                end	
            end  
                
            fhi = VecPhi[a,b]  
            eleminfalla = NeighbordCell[5] 
        
            vecreload = distmenorLLS_Asp(NeighbordCell,fhi)      
        
            contaBucleNorm = 0          
            for jk = b - 1 : b + 1 
                for ik = a - 1 : a + 1
                    contaBucleNorm += 1
                    VecPosi[ik,jk] = vecreload[contaBucleNorm]
                end	
            end  
        
            suma = 0
            sumarho = 0        
            VecPosiArray = VecPosi[2:Nboxx,2:Nboxy]
            B = filter(x -> x != -1.0, VecPosiArray)
            suma = sum(B)
            sumarho = sum(B.^rho)  
        
            tiempo[k] = 1.0/sumarho 
            acumt=acumt+tiempo[k]
            vectk1[k,:] = [k,acumt,tiempo[k],0,suma,sumarho,eleminfalla,a,b,fhi,iorigi,jorigi]
            contDisit = 0
            if VecPhi[a,b] != fhiFuera
                for ik = 1:length(vectFinSimul[:,1])
                    if a ==  vectFinSimul[ik,1] &&  b == vectFinSimul[ik,2]
                        break
                    else
                        contDisit +=1
                    end
                end
          
                if contDisit == length(vectFinSimul[:,1])
                    contFin +=1
                    if contFin > length(vectFinSimul[:,1])
                        break
                    else
                      vectFinSimul[contFin,1] = a
                      vectFinSimul[contFin,2] = b
                    end
                end
            end
        
            TempMean = 0  
            VecMeanTempo = zeros(Nboxx*Nboxy)
            for i = 1:length(VecPosiArray)
                if VecPosiArray[i]  >= 0.0
                    TempMean += 1
                    VecMeanTempo[TempMean] = VecPosiArray[i]
                end
            end
            
            vectparamestad[k,1] = mean(VecMeanTempo[1:TempMean])
            vectparamestad[k,2] = std(VecMeanTempo[1:TempMean])                
            VecMeanTempo[1:TempMean] = 0.0
            TempMean=0
            acumtref = acumt
            tiemporef = tiempo[k]
            indicenorm = indicenorm+1
            contNorm = contNorm+1
        end
      
              
    end  # close the WHILE loop
  
    MaxPasos = 0
    MaxTime = 0.0
    fin = 0
    MaxPasos=length(vectk1[:,1])-1
    fin=contNumAval    
    
    VecPosiFinal = VecPosi 
    return vectk1[1:MaxPasos,:],vectparamestad[1:length(vectparamestad[:,1]),:],VecPosiFinal
end
