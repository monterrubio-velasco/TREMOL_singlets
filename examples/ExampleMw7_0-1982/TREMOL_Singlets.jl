# Author:  Marisol Monterrubio-Velasco
# Contact: marisol.monterrubio@bsc.es

using PyPlot

YOURPATH = "userPATH"
include(joinpath("/home/mmonterr/Documentos/TREMOL/TREMOL_singlets/TREMOL_singlets/main/","FBM_Singlets.jl"))
include(joinpath("/home/mmonterr/Documentos/TREMOL/TREMOL_singlets/TREMOL_singlets/postprocessing/","postprocessing.jl"))
pathresults = "/home/mmonterr/Documentos/TREMOL/TREMOL_singlets/examples/ExampleMw7_0-1982/Results_PlotsFilesMw7_1982/"  

function TREMOL_Singlets(VecID,fhi_asp,fhi_bkg,strength_asp,DurTeo,VelAsp,Leff,Weff,VectorSeff,SaOri,YsizeOri,nbox,nk)   

# """
# TREMOL_Singlets.jl :  Main program that assigns the size and shape to the effective domain and asperity domain. Also this script gives the input values to the FBM algorithm. Lastly the output values goes to the postprocess function.
 
#Import function
# - include(joinpath("/YOURPATH/TREMOL_singlets/main/","FBM_Singlets.jl"))
# - include(joinpath("/YOURPATH/TREMOL_singlets/postprocessing/","postprocessing.jl"))

# Parameters
#  - `VecID::String`: Name to identify the simulated event.
#  - `fhi_asp::Float`: percentage of load to be transferred for any ruptured cell in the asperity domain
#  - `fhi_bkg::Float or Vector`: percentage of load to be transferred for any ruptured cell in the background
#  - `strength_asp::Integer`: strength value to define the asperity cells.
#  - `DurTeo::Float Vector`: rupture duration given by finite fault computation [seconds]
#  - `VelAsp :: Float Vector`: rupture velocity given by finite fault computation [km/s]
#  - `Leff::Float Vector`: Effective length size of the effective rupture area [km]
#  - `Weff::Float Vector`: Effective wide size of the effective rupture area [km]
#  - `VectorSeff::Float Vector`: Effective area [km²]
#  - `SaOri::Float Vector`: Ratio of the asperity size
#  - `YsizeOri::Float Vector`: Asperity area [km²]
#  - `nbox::Integer`: Number of cells assigned to each Seismic Source.
#  - `nk::Integer`: Number of times that will execute a same experiment to obtain statistics

# return
#   - `VecMagni::Vector`: Array that contains the results of the magnitude analysis coming from the function postprocessing.jl
# """
    
    rij = length(Leff)   # number of elements in the array 
    println("Effective Area -- Aeff=",Leff[rij].*Weff[rij],"km²--------------")
    IntervalAsp = 1
    lateralSizeGrid = nbox 
    totalNumberCells = lateralSizeGrid^2     # number of elements in the domain
    seffTotal = VectorSeff[rij]              # effective area size 
    nbox_x = round.(Int, sqrt(totalNumberCells/(Weff[rij]/Leff[rij])))    
    nbox_y = round.(Int,(Weff[rij]/Leff[rij])*nbox_x)    
    println("Grid-----Nx =",nbox_x,"-----Ny=",nbox_y,"-----")
    cellSize = seffTotal/totalNumberCells       # Individual cell size
    println(cellSize)    
    smin = 0    
#   vector to save the magnitude results   
    VecMagni = zeros(nk+1,24)   
#   loop to realize statistics using same input data but different asperity size
    for ik = 1:nk              
        srand(ik+1)
        smin = round(Int,totalNumberCells*1.5)    
        a = rand()     
        VectorSaOri = SaOri[rij] + (a*(SaOri[rij]/2))
        println("VectorSaOri",VectorSaOri)
        smin += round.(Int,(a*smin/2))
               
        #################################################################### 
#       Number of cells in each asperities
        VectorCellsSaOri = round.(Int,totalNumberCells*VectorSaOri)  
#       Asperity lateral cell size if is considered as square                
        VectorSaOriLateralSize = round.(Int,sqrt.(VectorCellsSaOri))  
#       Asperity lateral cell size at Y-axis        
        VectorSaOriLateralSize_y = round(Int,sqrt(VectorCellsSaOri./(Weff[rij]./Leff[rij]))) 
#       Asperity lateral cell size at X-axis       
        VectorSaOriLateralSize_x = round(Int,(Weff[rij]/Leff[rij])*VectorSaOriLateralSize_y)          
#       Each seismic source in the array has same nbox_y so we compute their nbox_x
        VectorCellSeffArray = round.(Int,totalNumberCells./nbox_y)
#       Each Asperity is located in the center of each square
        VectorCenterAperities = zeros(1,2)
        cumY = 0
        for i = 1:length(VectorCenterAperities[:,1])     
            VectorCenterAperities[i,1] = VectorCellSeffArray[i]/2 + cumY
            VectorCenterAperities[i,2] = nbox_y/2
            println(VectorCenterAperities[i,:])  
            cumY += VectorCellSeffArray[i] 
        end
#       Coordinates that defines each asperirties area
        VectorCoordsAperities_x = zeros(1,4)
        VectorCoordsAperities_y = zeros(1,4)
    
        VectorCoordsAperities_x[1,1] = (VectorCenterAperities[1,1]-VectorSaOriLateralSize_x[1]/2)
        VectorCoordsAperities_y[1,1] = (VectorCenterAperities[1,2]+VectorSaOriLateralSize_y[1]/2)              
        VectorCoordsAperities_x[1,2] = (VectorCenterAperities[1,1]-VectorSaOriLateralSize_x[1]/2)
        VectorCoordsAperities_y[1,2] = (VectorCenterAperities[1,2]-VectorSaOriLateralSize_y[1]/2)          
        VectorCoordsAperities_x[1,3] = (VectorCenterAperities[1,1]+VectorSaOriLateralSize_x[1]/2)
        VectorCoordsAperities_y[1,3] = (VectorCenterAperities[1,2]-VectorSaOriLateralSize_y[1]/2)              
        VectorCoordsAperities_x[1,4] = (VectorCenterAperities[1,1]+VectorSaOriLateralSize_x[1]/2)
        VectorCoordsAperities_y[1,4] = (VectorCenterAperities[1,2]+VectorSaOriLateralSize_y[1]/2)    

        VectorCoordsAperities_x = round.(Int,VectorCoordsAperities_x)
        VectorCoordsAperities_y = round.(Int,VectorCoordsAperities_y)
     
#       Construction of the matrix 
        MatrizLoads = zeros(nbox_x,nbox_y)
        MatrizStrengt = ones(nbox_x,nbox_y)
        MatrizPhi = zeros(nbox_x,nbox_y) 
        for j = 1: nbox_y 
            for i = 1: nbox_x    
                MatrizLoads[i,j] = rand()        
                MatrizPhi[i,j] = fhi_bkg                
                if (VectorCenterAperities[1,1]-VectorSaOriLateralSize_y[1]/2) <= i <= (VectorCenterAperities[1,1]+ VectorSaOriLateralSize_y[1]/2) && (VectorCenterAperities[1,2]-VectorSaOriLateralSize_x[1]/2) <= j <= (VectorCenterAperities[1,2]+VectorSaOriLateralSize_x[1]/2) 
                    MatrizStrengt[i,j] = rand(strength_asp-IntervalAsp:strength_asp+IntervalAsp)
                    MatrizPhi[i,j] = fhi_asp
                end                          
            end
        end
          
        maxiStrength = maximum(MatrizStrengt)
        minStrength = minimum(MatrizStrengt)    
        MatrizStrengtInitial = copy(MatrizStrengt)

#       Plotting the strength domains of the asperity and the background
        if ik == 1
            PyPlot.figure(1,(10,15))
            img = PyPlot.imshow(MatrizStrengt,cmap="BrBG")
            cbar=PyPlot.colorbar(img,ticks=[0.0:0.5:maxiStrength])         
            PyPlot.title(VecID[rij]*L"$-S_{syn}$= "*string(round(VectorSaOri,2)),fontsize = 15)       
            PyPlot.xlabel("cells N",fontsize = 17,weight="semibold")
            PyPlot.ylabel("cells N",fontsize = 17,weight="semibold")
            PyPlot.xticks(fontsize=15,weight="semibold")
            PyPlot.yticks(fontsize=15,weight="semibold")  
            PyPlot.text((nbox_x+20), round(Int,(nbox_y/2)+10), L"$\gamma(i,j)$", color="k",fontsize = 8,rotation=90, va="bottom", ha="center")
            PyPlot.savefig(pathresults*VecID[rij]*"AsperezaSpatial.pdf",dpi=400)
            PyPlot.hold(false)   
        end
              
#       Program that computes FBM algorithm
        NumberAsperities = 1
          
        datos,vecEstadistico,VecPosiFinal= FBM_Singlets(nbox_x,nbox_y,smin,MatrizLoads,MatrizStrengtInitial,MatrizStrengt,fhi_bkg,fhi_asp,MatrizPhi)
        
#       Program that realize the post-process analysis               
        VecMagniLoop = postprocessing(ik,datos,VectorCenterAperities,VectorSaOriLateralSize_y,VectorSaOriLateralSize_x,VectorCoordsAperities_x,VectorCoordsAperities_y,nbox_x,nbox_y, smin,vecEstadistico,VectorSaOri,SaOri[rij],totalNumberCells,a,cellSize,VelAsp[rij],Weff[rij],Leff[rij],DurTeo[rij],pathresults,VecID[rij])
              
#       Results of the post processing save in the results vector
        VecMagni[ik,:] = VecMagniLoop           
              
    end
    
    VecMagni[nk+1,1] = strength_asp
    VecMagni[nk+1,4] = fhi_asp
    VecMagni[nk+1,2] = fhi_bkg
    VecMagni[nk+1,3] = lateralSizeGrid
    VecMagni[nk+1,5] = cellSize   
  
#   Save the results vector
    writedlm(pathresults*VecID[rij]*"-MagnitudeStatisticalResults.dat",VecMagni)

end
