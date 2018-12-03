# Author:  Marisol Monterrubio-Velasco
# Contact: marisol.monterrubio@bsc.es
# GNU GENERAL PUBLIC LICENSE
# ==========================
# Version 3, 29 June 2007
# Copyright © 2007 Free Software Foundation, Inc. <http://fsf.org/>

include("gutenberRichAspDef.jl")
include("bmemag.jl")

function calcuMagniSpaceTimeSinglets(ik,Y,Nbox_x,Nbox_y,smin,VectorCenterAperities,VectorSaOriLateralSize_y,VectorSaOriLateralSize_x, CellSize,velAsp, Weff,Leff,pathresults,VecID)

# """
# calcuMagniSpaceTimeMultiSinglets.jl : Main program that assigns the size and shape to the effective domain and asperity domain. Also this script gives the input values to the FBM algorithm. Lastly the output values goes to the postprocess function.
# 
#Import function
# - include(joinpath("/YOURPATH/TREMOL_singlets/postprocessing/","gutenberRichAspDef.jl"))
# - include(joinpath("/YOURPATH/TREMOL_singlets/postprocessing/","bmemag.jl"))

#Parameters
# - `Y::Array`:raw data coming from the FBM algorithm. This data base contains the rupture information of model  
## Each row contains: 
    ## 1. k:number of step
    ## 2. acumt: tiempo acumulado (T_k = sum(tiempo[1:k]))
    ## 3. tiempo[k]: inter-event time [dimensionless]
    ## 4. (0 or 1): identifier to normal or avalanche event
    ## 5. suma: sum of the load in all the cells
    ## 6. sumarho= 1/tiempo[k]  (rupture rate)
    ## 7. parametrosigma: load value of the cell chosen to fail
    ## 8. a: coordinate in the X-axis, of the cell chosen to fail
    ## 9. b: coordinate in the Y-axis, of the cell chosen to fail
    ##10. fhi: load-transfer value
    ##11-12. (iorigi,jorigi): coordinates of the cell chosen to fail at the first stage of the searching algortithm defined in the function contaravalan_Asp. 
    
# - `Nbox_x::Integer`: number of cells in X-axis of the domain \Omega
# - `Nbox_y::Integer`: number of cells in Y-axis of the domain \Omega`
# - `smin::Integer`: number of steps that realize the algorithm 
# - `VectorCenterAperities::Vector Integer`: size(1), Central coordinates of the aspertity (x,y)
# - `VectorSaOriLateralSize_y::Vector`: Asperity lateral size in the Y-axis
# - `VectorSaOriLateralSize_x::Vector`: Asperity lateral size in the X-axis
# - `VectorCoordsAperities_x::Vector Integer`: coordinates of the asperity vertex in the X-axis
# - `VectorCoordsAperities_y::Vector Integer`: coordinates of the asperity vertex in the Y-axis 
# - `CellSize::Float`: Size of a cell in km²`
# - `velAsp::Float`: rupture velocity
# - `Weff::Float`: width of the effective area in km 
# - `Leff::Float`: length of the effective area in km 

#  return
# - `vecMagni::Array`: matrix that contains the results of the magnitude analysis
# - VecNewAvalSpaceTime::Array`: matrix that contains the data of the regrouping algorithm coming from the Y matrix  
## Each row contains: 
# - VecNewAvalSpaceTime[contNumAvalSpaceTime,1] = contadorAval: number of elements considering in the  
#                                                   regrouping algorithm
# - VecNewAvalSpaceTime[contNumAvalSpaceTime,2] = original number in the raw catalog Y
# - VecNewAvalSpaceTime[contNumAvalSpaceTime,3] = cumulative time (T_k = sum(tiempo[1:k])) computed in 
#                                                   the raw catalog Y
# - VecNewAvalSpaceTime[contNumAvalSpaceTime,4] = inter-event time [dimensionless] computed in the raw
#                                                   catalog Y
# - VecNewAvalSpaceTime[contNumAvalSpaceTime,5] = normal-avalanche code considered in the new regrouping
# - VecNewAvalSpaceTime[contNumAvalSpaceTime,6] = x coordinate
# - VecNewAvalSpaceTime[contNumAvalSpaceTime,7] = y coordinate
# - VecNewAvalSpaceTime[contNumAvalSpaceTime,8] = rupture rate  
# """

    AreaSUB = CellSize 
    DurSingCell = sqrt(CellSize)/velAsp
    VecNew = copy(Y)
    contador = 0
    VecNewAval = zeros(length(VecNew[:,1]),6)
    vecMagni = zeros(4,16)

# Avalanche events are separated from the data raw 
    for i=1:length(VecNew[:,1]) 
        if VecNew[i,4] == 1
            contador += 1     
            VecNewAval[contador,1] = VecNew[i,1]  #cumulative number of event
            VecNewAval[contador,2] = VecNew[i,2]  #cumulative time
            VecNewAval[contador,3] = VecNew[i,3]  #inter event time
            VecNewAval[contador,4] = VecNew[i,4]  #avalanche or normal identifier
            VecNewAval[contador,5] = VecNew[i,8]  #X-position
            VecNewAval[contador,6] = VecNew[i,9]  #Y-position        
        end     
    end  
  
    rdist = 0.0
    rmin = sqrt(2) # the avalanches are cluster considering the spatial neareast neighbors   
    VecNewAvalSpaceTime = zeros(contador,8)
    VecAvalanchas = zeros(contador)
    VecTimeAvalInicial = zeros(contador)
    VecTimeAvalFinal = zeros(contador)
    VecTimeRuptureAvalanche = zeros(contador)
    contadorAval = 0
    contNumAvalSpaceTime = 1
    contAva = 0
  
    for j=2:contador  
        rdist = sqrt((VecNewAval[j,5]-VecNewAval[j-1,5])^2 + (VecNewAval[j,6]-VecNewAval[j-1,6])^2)
               
        if (rdist <= rmin) || (VectorCenterAperities[1,1]-VectorSaOriLateralSize_y[1]/2) <= VecNewAval[j,5] <= (VectorCenterAperities[1,1]+ VectorSaOriLateralSize_y[1]/2) && (VectorCenterAperities[1,1]-VectorSaOriLateralSize_y[1]/2) <= VecNewAval[j-1,5] <= (VectorCenterAperities[1,1]+ VectorSaOriLateralSize_y[1]/2) || (VectorCenterAperities[1,2]-VectorSaOriLateralSize_x[1]/2) <= VecNewAval[j,6] <= (VectorCenterAperities[1,2]+ VectorSaOriLateralSize_x[1]/2) && (VectorCenterAperities[1,2]-VectorSaOriLateralSize_x[1]/2) <= VecNewAval[j-1,6] <= (VectorCenterAperities[1,2]+ VectorSaOriLateralSize_x[1]/2) 

            contadorAval += 1   
            contNumAvalSpaceTime += 1
            VecNewAvalSpaceTime[contNumAvalSpaceTime,1] = contadorAval
            VecNewAvalSpaceTime[contNumAvalSpaceTime,2] = VecNewAval[j,1] 
            VecNewAvalSpaceTime[contNumAvalSpaceTime,3] = VecNewAval[j,2]
            VecNewAvalSpaceTime[contNumAvalSpaceTime,4] = VecNewAval[j,3]  # inter-event time
            VecNewAvalSpaceTime[contNumAvalSpaceTime,5] = VecNewAval[j,4]
            VecNewAvalSpaceTime[contNumAvalSpaceTime,6] = VecNewAval[j,5]   # x-position
            VecNewAvalSpaceTime[contNumAvalSpaceTime,7] = VecNewAval[j,6]   # y-position   
            VecNewAvalSpaceTime[contNumAvalSpaceTime,8] = 1/VecNewAval[j,3] 
            if contadorAval == 1
                contAva += 1
                VecTimeAvalInicial[contAva] = VecNewAvalSpaceTime[contNumAvalSpaceTime,3]
            end
            if j == contador
                VecAvalanchas[contAva] = contadorAval
            end
      
        else                     
            contadorAval = 0         
            contNumAvalSpaceTime += 1
            VecNewAvalSpaceTime[contNumAvalSpaceTime,1] = contadorAval
            VecNewAvalSpaceTime[contNumAvalSpaceTime,2] = VecNewAval[j,1] 
            VecNewAvalSpaceTime[contNumAvalSpaceTime,3] = VecNewAval[j,2]
            VecNewAvalSpaceTime[contNumAvalSpaceTime,4] = VecNewAval[j,3]
            VecNewAvalSpaceTime[contNumAvalSpaceTime,5] = 0
            VecNewAvalSpaceTime[contNumAvalSpaceTime,6] = VecNewAval[j,5] # x
            VecNewAvalSpaceTime[contNumAvalSpaceTime,7] = VecNewAval[j,6] # y
            VecNewAvalSpaceTime[contNumAvalSpaceTime,8] = 1/VecNewAval[j,3]       
            if VecNewAvalSpaceTime[contNumAvalSpaceTime-1,1] != 0
                VecAvalanchas[contAva] = VecNewAvalSpaceTime[contNumAvalSpaceTime-1,1] 
                VecTimeAvalFinal[contAva] = VecNewAvalSpaceTime[contNumAvalSpaceTime,3]
                VecTimeRuptureAvalanche[contAva] =  VecTimeAvalFinal[contAva] - VecTimeAvalInicial[contAva]  
            end
            if j == contador
                VecAvalanchas[contAva] = VecNewAvalSpaceTime[contNumAvalSpaceTime-1,1] 
                VecTimeAvalFinal[contAva] = VecNewAvalSpaceTime[contNumAvalSpaceTime,3]
            end
        end     
    end
    
    r=0
    VecAval = zeros(contador)
    VecArea = zeros(contador)
    VecTimeIni = zeros(contador)
    VecTimeFin = zeros(contador)
    DuraciAval = zeros(contador)
    VecRateRuptAval = zeros(contador)
    maxAv = 0
    maxTiIni = 0
    maxTiFin = 0
    maxArea = 0
    maxDurat = 0
    
    for i=1:contAva
        if VecAvalanchas[i] >= 2
            r += 1      
            VecAval[r] = VecAvalanchas[i]
            VecArea[r] = AreaSUB * VecAvalanchas[i]
            VecTimeIni[r] = VecTimeAvalInicial[i]
            VecTimeFin[r] = VecTimeAvalFinal[i]
            DuraciAval[r] = sqrt(VecArea[r])/velAsp  
            VecRateRuptAval[r] = VecAval[r] / VecTimeRuptureAvalanche[r]  # rupture velocity (dimensionless)      
            if VecAval[r] > maxAv
                maxAv = VecAval[r]
                maxTiIni = VecTimeIni[r]
                maxTiFin = VecTimeFin[r]
                maxArea = VecArea[r]
                maxDurat = DuraciAval[r] 
            end
        end
    end
  
    VecVelocidad =   VecNewAvalSpaceTime[1:contNumAvalSpaceTime,8]
    maxiVelocidad = maximum(VecVelocidad)
    VectorSaLateralSizeAsp_x = (Weff/Leff)*sqrt(maxArea)
    VectorSaLateralSizeAsp_y = maxArea/VectorSaLateralSizeAsp_x
    maxDuratRect = VectorSaLateralSizeAsp_y/velAsp
    maxDuratSq = sqrt(maxArea)/velAsp

# maximum rupture percentage 
    maxAval = maximum(VecAval[1:r])
    B = sort(VecAval[1:r])
    maxiArea = maximum(VecArea[1:r])
    
# scale relations to compute magnitud (Q. Rodriguez and Ottemöller, 2013) 

  # ==============..............Somerville..............=============
    
    VecMagSomerville = (log10.(VecArea[1:r]) + 4.393 ) / 0.991 
    MaxMagnitud = (log10.(maxArea) + 4.393 ) / 0.991 
    maxMag = maximum(VecMagSomerville)	
 
    paramGRSomer, cuenrepSomer, vecCumMagniSomer = gutenberRichAspDef(ik,VecMagSomerville,"Somer",AreaSUB,pathresults*VecID)
    meanmSomer, bSomer, sigSomer, av2Somer = bmemag(VecMagSomerville)
    NumberValueBSomer = cuenrepSomer/length(VecMagSomerville)	
  
    vecMagni[1,1] = meanmSomer
    vecMagni[1,2] = bSomer
    vecMagni[1,3] = sigSomer
    vecMagni[1,4] = maximum(VecMagSomerville)
    vecMagni[1,5] = NumberValueBSomer  
    vecMagni[1,6] = maxAval # maximum Avalanche
    vecMagni[1,7] = MaxMagnitud # maximum Magnitude
    vecMagni[1,8] = maxArea/velAsp #tiempo de ruptura [segundos]
    vecMagni[1,9] = maxDuratRect
    vecMagni[1,10] = maxDuratSq #qSomerTelesca
    vecMagni[1,11] = maxiArea #minSomerTelesca
    vecMagni[1,12] = maxDurat
    vecMagni[1,13] = cuenrepSomer
    vecMagni[1,14] = length(VecMagSomerville)
    vecMagni[1,15] = av2Somer	
  
  # ==============..............Mai..............=============
    VecMai_Large = (log10.(VecArea[1:r]) + 5.581 ) / 1.137	
    paramMaiLS, cuenrepMaiLS, vecCumMagniMaiLS = gutenberRichAspDef(ik,VecMai_Large,"Mai-L",AreaSUB,pathresults*VecID)
    meanmMaiLS, bMaiLS, sigMaiLS, av2MaiLS = bmemag(VecMai_Large)
    NumberValueBMaiLS = cuenrepMaiLS/length(VecMai_Large)
  
    vecMagni[2,1] = meanmMaiLS
    vecMagni[2,2] = bMaiLS
    vecMagni[2,3] = sigMaiLS
    vecMagni[2,4] = maximum(VecMai_Large)
    vecMagni[2,5] = NumberValueBMaiLS
    vecMagni[2,6] = 0 
    vecMagni[2,7] = 0 
    vecMagni[2,8] = 0 
    vecMagni[2,9] = 0 
    vecMagni[2,10] = 0 
    vecMagni[2,11] = 0 
    vecMagni[2,12] = 0  
    vecMagni[2,13] = cuenrepMaiLS
    vecMagni[2,14] = length(VecMai_Large)
    vecMagni[2,15] = av2MaiLS

# ==============..............Mai-VeryLarge Asperity criteria..............=============
    VecMai_VeryLarge = (log10.(VecArea[1:r]) + 6.013 ) / 1.146 	
    paramMaiVLS, cuenrepMaiVLS, vecCumMagniMaiVLS = gutenberRichAspDef(ik,VecMai_VeryLarge,"Mai-VL",AreaSUB,pathresults*VecID)
    meanmMaiVLS, bMaiVLS, sigMaiVLS, av2MaiVLS = bmemag(VecMai_VeryLarge)
    NumberValueBMaiVLS = cuenrepMaiVLS/length(VecMai_VeryLarge) 	
  
    vecMagni[3,1] = meanmMaiVLS
    vecMagni[3,2] = bMaiVLS
    vecMagni[3,3] = sigMaiVLS
    vecMagni[3,4] = maximum(VecMai_VeryLarge)
    vecMagni[3,5] = NumberValueBMaiVLS
    vecMagni[3,6] = 0
    vecMagni[3,7] = 0
    vecMagni[3,8] = 0 
    vecMagni[3,9] = 0
    vecMagni[3,10] = 0 
    vecMagni[3,11] = 0 
    vecMagni[3,12] = 0  
    vecMagni[3,13] = cuenrepMaiVLS
    vecMagni[3,14] = length(VecMai_VeryLarge)
    vecMagni[3,15] = av2MaiVLS

# ==============..............Ramirez..............====================================
    VecRamirez = (2/3 * log10.( (VecArea[1:r]/(7.78*1.0e-9)).^(1/0.550)) ) - 6.07
    paramRamiAsper, cuenrepRamiAsper, vecCumMagniRamiAsper = gutenberRichAspDef(ik,VecRamirez,"Ramirez",AreaSUB,pathresults*VecID)
    meanmRami, bRami, sigRami, av2Rami = bmemag(VecRamirez)
    NumberValueBRamiAsper = cuenrepRamiAsper/length(VecRamirez)
	
    vecMagni[4,1] = meanmRami
    vecMagni[4,2] = bRami
    vecMagni[4,3] = sigRami
    vecMagni[4,4] = maximum(VecRamirez)
    vecMagni[4,5] = NumberValueBRamiAsper
    vecMagni[4,6] = 0 
    vecMagni[4,7] = 0 
    vecMagni[4,8] = 0 
    vecMagni[4,9] = 0 
    vecMagni[4,10] = 0 
    vecMagni[4,11] = 0 
    vecMagni[4,12] = 0  
    vecMagni[4,13] = cuenrepMaiVLS
    vecMagni[4,14] = length(VecRamirez)
    vecMagni[4,15] = av2Rami 
    
    
#   Plots 

#   plot Magnitude vs Time of the simulated earthquakes   
    PyPlot.figure(51,(15,10)) 
    PyPlot.hold("True")
    PyPlot.semilogx(VecTimeIni[1:r],VecMagSomerville[1:r],"o:",lw = 0.5,ms=2,label="Somerville scale relation")
    PyPlot.semilogx(maxTiIni, MaxMagnitud, "ro", ms=4)
    PyPlot.xticks(fontsize=15,weight="semibold")
    PyPlot.yticks(fontsize=15,weight="semibold")  
    xlabel("Time [dimensionless]", fontsize = 17,weight="semibold")
    ylabel("Magnitude",fontsize = 17,weight="semibold")
    if ik == 1
      PyPlot.legend(fontsize=14, markerscale=2)
    end
    PyPlot.savefig(pathresults*VecID*"-MagnitudeTime.pdf",dpi=200)
    PyPlot.hold("False")
# ==================================================================================================
    
#   Magnitude Histogram     
    PyPlot.figure(17,(10,10))
    h = plt[:hist](VecMagSomerville[1:r],22,label="Somerville_scale_relation")
    xlim(4,9)
    PyPlot.xticks(fontsize=15,weight="semibold")
    PyPlot.yticks(fontsize=15,weight="semibold")  
    ylabel("Frequency",fontsize = 17,weight="semibold")
    xlabel("Magnitude",fontsize = 17,weight="semibold")
    if ik == 1
      PyPlot.legend(fontsize=14, markerscale=2)
    end
    PyPlot.savefig(pathresults*VecID*"-FrequencyMagnitude.pdf",dpi=200)
# ==================================================================================================
      
#   Rupture Duration vs Number of simulated earthquake 
    PyPlot.figure(59,(10,10))
    plot(collect(1:r),DuraciAval[1:r],"--o",label="Rupture_Duration_of_Earthquakes")
    PyPlot.xticks(fontsize=15,weight="semibold")
    PyPlot.yticks(fontsize=15,weight="semibold")  
    xlabel("Number of Avalanches",fontsize = 17,weight="semibold")
    ylabel("Rupture Duration [s]",fontsize = 17,weight="semibold")
    if ik == 1
      PyPlot.legend(fontsize=14, markerscale=2)
    end
    PyPlot.savefig(pathresults*VecID*"-DurationFrequency.pdf",dpi=200)
 # ==================================================================================== 
   
    return vecMagni,VecNewAvalSpaceTime[1:contNumAvalSpaceTime,:]

end
