# Author:  Marisol Monterrubio-Velasco
# Contact: marisol.monterrubio@bsc.es
# GNU GENERAL PUBLIC LICENSE
# ==========================
# Version 3, 29 June 2007
# Copyright © 2007 Free Software Foundation, Inc. <http://fsf.org/>

include("calcuMagniSpaceTimeSinglets.jl")
include("plotcoordenadasSingletes.jl")
   
function postprocessing(ik,datos,VectorCenterAperities,VectorSaOriLateralSize_y ,VectorSaOriLateralSize_x,VectorCoordsAperities_x,VectorCoordsAperities_y,nbox_x,nbox_y, smin,vecEstadistico,VectorSaOri,SaOri,TotalNumberCells,a,CellSize,velAsp,Weff,Leff,DurTeo,pathresults,VecID)

# """
# postprocessing.jl : Main program that assigns the size and shape to the effective domain and asperity domain. Also this script gives the input values to the FBM algorithm. Lastly the output values goes to the postprocess function.
# 
#Import function
# - include(joinpath("/YOURPATH/TREMOL_singlets/main/","calcuMagniSpaceTimeMultiSinglets.jl"))
# - include(joinpath("/YOURPATH/TREMOL_singlets/postprocessing/","plotcoordenadasSingletesbis.jl"))

#Parameters
# - `datos::Array`: size(smin,12) raw data coming from the FBM algorithm. This data base contains the rupture information of model  
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
    
# - `VectorCenterAperities::Vector Integer`: size(1), Central coordinates of the aspertity (x,y)
# - `VectorSaOriLateralSize_y::Vector`: Asperity lateral size in the Y-axis
# - `VectorSaOriLateralSize_x::Vector`: Asperity lateral size in the X-axis
# - `VectorCoordsAperities_x::Vector Integer`: coordinates of the asperity vertex in the X-axis
# - `VectorCoordsAperities_y::Vector Integer`: coordinates of the asperity vertex in the Y-axis 
# - `nbox_x::Integer`: number of cells in X-axis of the domain  \Omega
# - `nbox_y::Integer`: number of cells in Y-axis of the domain \Omega 
# - `smin::Integer`: number of steps that realize the algorithm 
# - `vecEstadistico::Array`:  size(smin+10,2). Mean and Standard deviation of the load in the system, computed at each time step considering only the cells active
# - `VectorSaOri::Float Vector`: Random size  of the asperity 
# - `SaOri::Float Vector`: Original Ratio of the asperity size
# - `TotalNumberCells::Integer`: number of cells in the domain \Omega (nbox_x * nbox_y)
# - `a::Float`: random number
# - `CellSize::Float`: Size of a cell in km²
# - `velAsp::Float`: rupture velocity
# - `Weff::Float`: width of the effective area in km 
# - `Leff::Float`: length of the effective area in km 
# - `DurTeo::Float`: rupture duration


#return
# - `vecMagniLoop::Array`: results of tha magnitude analysis 
## Each row contains: 
    ## 1. vecMagniLoop[1,1]: b-value computed trough the function bmemag.jl using Somerville relation
    ## 2. vecMagniLoop[1,2]: maximum magnitude using Somerville relation
    ## 3. vecMagniLoop[1,3]: b-value computed trough the function bmemag.jl using Mai relation
    ## 4. vecMagniLoop[1,4]: maximum magnitude using Mai relation
    ## 5. vecMagniLoop[1,5]: b-value computed trough the function bmemag.jl using Mai-VL relation
    ## 6. vecMagniLoop[1,6]: maximum magnitude using Mai-VL relation
    ## 7. vecMagniLoop[1,7]: b-value computed trough the function bmemag.jl using Ramirez relation
    ## 8. vecMagniLoop[1,8]: maximum magnitude using Ramirez relation
    ## 9. vecMagniLoop[1,9]: ratio of the largest simulated earthquake [cells] and the total number of cells
    ##10. vecMagniLoop[1,10]: largest simulated earthquake in [cells]
    ##11. vecMagniLoop[1,11]: VectorSaOri
    ##12. vecMagniLoop[1,12]: a
    ##13. vecMagniLoop[1,13]: SaOri
    ##14. vecMagniLoop[1,14]: smin
    ##15. vecMagniLoop[1,15]: maximum magnitude using Somerville relation
    ##16. vecMagniLoop[1,16]: length of the area earhquake (in km^2) divided by the rupture velocity (km/s). 
    ##17. vecMagniLoop[1,17]: equivalent rupture time in seconds, considering the longest length divided by the rupture velocity
    ##18. vecMagniLoop[1,18]: equivalent rupture time in seconds considering the root square of the area divided by the rupture velocity
    ##19. vecMagniLoop[1,19]: rupture velocity, velAsp
    ##20. vecMagniLoop[1,20]: area of the largest simulated earthquake
    ##21. vecMagniLoop[1,21]: sqrt(area of the largest simulated earthquake)/velAsp
    ##22. vecMagniLoop[1,22]: sqrt(area of the largest simulated earthquake)/DurTeo
    ##23. vecMagniLoop[1,23]: velAsp*DurTeo
    ##24. vecMagniLoop[1,24]: sqrt(area of the largest simulated earthquake)/velAsp 
# """
 
    vecMagniAval,vecNewAvalSpaceTime = calcuMagniSpaceTimeSinglets(ik,datos,nbox_x,nbox_y,smin,VectorCenterAperities,VectorSaOriLateralSize_y, VectorSaOriLateralSize_x, CellSize,velAsp, Weff,Leff,pathresults,VecID)
     
    plotcoordenadasSingletes(ik,vecNewAvalSpaceTime,length(datos[:,1]),nbox_x,nbox_y,VectorCenterAperities,VectorSaOriLateralSize_y,VectorSaOriLateralSize_x,VectorCoordsAperities_x,VectorCoordsAperities_y, CellSize,vecEstadistico,datos,pathresults,VecID)
      
    vecMagniLoop = zeros(1, 24)
    vecMagniLoop[1, 1] = vecMagniAval[1,2]
    vecMagniLoop[1, 2] = vecMagniAval[1,4] 
    vecMagniLoop[1, 3] = vecMagniAval[2,2]
    vecMagniLoop[1, 4] = vecMagniAval[2,4] 
    vecMagniLoop[1, 5] = vecMagniAval[3,2]
    vecMagniLoop[1, 6] = vecMagniAval[3,4] 
    vecMagniLoop[1, 7] = vecMagniAval[4,2]
    vecMagniLoop[1, 8] = vecMagniAval[4,4]    
    vecMagniLoop[1, 9] = vecMagniAval[1,6]/TotalNumberCells 
    vecMagniLoop[1, 10] = vecMagniAval[1,6]
    vecMagniLoop[1, 11] = VectorSaOri
    vecMagniLoop[1, 12] = a
    vecMagniLoop[1, 13] = SaOri
    vecMagniLoop[1, 14] = smin
    vecMagniLoop[1, 15] = vecMagniAval[1,7]
    vecMagniLoop[1, 16] = vecMagniAval[1,8]
    vecMagniLoop[1, 17] = vecMagniAval[1,9]
    vecMagniLoop[1, 18] = vecMagniAval[1,10]
    vecMagniLoop[1, 19] = velAsp
    vecMagniLoop[1, 20] = vecMagniAval[1,11]
    vecMagniLoop[1, 21] = sqrt(vecMagniAval[1,11])/velAsp
    vecMagniLoop[1, 22] = sqrt(vecMagniAval[1,11])/DurTeo
    vecMagniLoop[1, 23] = velAsp*DurTeo
    vecMagniLoop[1, 24] = vecMagniAval[1,12]
              
       
    return vecMagniLoop               
    
end
