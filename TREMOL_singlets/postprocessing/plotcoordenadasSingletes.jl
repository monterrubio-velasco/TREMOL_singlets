# Author:  Marisol Monterrubio-Velasco
# Contact: marisol.monterrubio@bsc.es
# GNU GENERAL PUBLIC LICENSE
# ==========================
# Version 3, 29 June 2007
# Copyright © 2007 Free Software Foundation, Inc. <http://fsf.org/>

function plotcoordenadasSingletes(ik,datos,smin,Nboxx,Nboxy,VectorCenterAperities,VectorSaLateralSize_y,VectorSaLateralSize_x,VectorCoordsAperities_x,VectorCoordsAperities_y, CellSize,vecEstadistico,datosOri,pathresults,VecID)

# """
# plotcoordenadasSingletes.jl. This function plot the results of spatial distributiobn, inter-event rate, mean load, 
# 

#Parameters
# - `datos::Array`: Vector that contains the regrouping data coming from calcuMagniSpaceTimeMultiSinglets.jl function. This vector is in fact the data of simulated earthaquakes
# - `smin::Integer`: number of steps that realize the algorithm 
# - `Nboxx::Integer`: number of cells in X-axis of the domain  \Omega
# - `Nboxy::Integer`: number of cells in Y-axis of the domain \Omega 
# - `VectorCenterAperities::Vector Integer`: size(1), Central coordinates of the aspertity (x,y)
# - `VectorSaOriLateralSize_y::Vector`: Asperity lateral size in the Y-axis
# - `VectorSaOriLateralSize_x::Vector`: Asperity lateral size in the X-axis
# - `VectorCoordsAperities_x::Vector Integer`: coordinates of the asperity vertex in the X-axis
# - `VectorCoordsAperities_y::Vector Integer`: coordinates of the asperity vertex in the Y-axis 
# - `CellSize::Float`: Size of a cell in km²
# - `vecEstadistico::Array`:  size(smin+10,2). Mean and Standard deviation of the load in the system, computed at each time step considering only the cells active
# -`datos::Array`: raw data coming from the FBM algorithm. This data base contains the rupture information of model 
# """

#   Plot the mean load value vs time 
    PyPlot.figure(11234,(10,10))
    PyPlot.plot(log10.(datosOri[1:length(vecEstadistico[:,1]),2]),vecEstadistico[1:length(vecEstadistico[:,1]),1],"-b",label="Mean_Load")
    PyPlot.xlabel("Time [dimensionless]",fontsize = 17,weight="semibold")
    PyPlot.ylabel("Mean Load",fontsize = 17,weight="semibold")
    PyPlot.xticks(fontsize=15,weight="semibold")
    PyPlot.yticks(fontsize=15,weight="semibold") 
    if ik == 1
      PyPlot.legend(fontsize=14, markerscale=2)
    end
    PyPlot.savefig(pathresults*VecID*"-Mean_Load.pdf",dpi=200)
#   =======================================================================================  

    Z=zeros(1,length(datos[1,:]))
    datos = [datos;Z]
    smin=length(datos[:,1])
    x1=datos[2:smin,6] 
    y1=datos[2:smin,7]
    avalNorm=datos[2:smin,5]
    XX=zeros(length(x1)+1,3)
    l1=zeros(length(x1)+1)
    TimeSerie=datos[2:smin,3]
    frontx=Nboxx+20
    fronty=Nboxy+20
    VecSpatialAva=zeros(length(x1)+1,3)
    
    srand(20)
    xmin=0
    xmax=frontx
    ymin=0
    ymax=fronty
    v = [xmin,xmax,ymin,ymax]
    firstAval=0
    XfirstAval=0
    YfirstAval=0
    for i=1:length(x1)-1
        if avalNorm[i]==1
            firstAval+=1
            if firstAval==1
                XfirstAval=x1[i]
                YfirstAval=y1[i]
            end
        end
    end


#  Inter-event rate of Avalanche (in red color) and Normal events (in blue color) 
    XsecondAval=0
    YsecondAval=0
    AvalNumber=0
    Velocidad=zeros(firstAval,2)
    TiempoAval=zeros(firstAval,3)
    TimeAval=0
    xx1=0
    yy1=0
    for i=1:(length(x1)-2)
        if avalNorm[i]==1
            AvalNumber+=1
            xx1=x1[i]
            yy1=y1[i]    
            if (VectorCenterAperities[1,1]-VectorSaLateralSize_y[1]/2) <= xx1 <= (VectorCenterAperities[1,1]+ VectorSaLateralSize_y[1]/2) && (VectorCenterAperities[1,2]-VectorSaLateralSize_x[1]/2) <= yy1 <= (VectorCenterAperities[1,2]+ VectorSaLateralSize_x[1]/2) 
                nflag=1
            else     
                nflag=0
            end
            if  x1[i] == XfirstAval && y1[i] == YfirstAval
                DistanciaAval = 0
                Velocidad[AvalNumber,1]=0
                TiempoAval[AvalNumber,1]= TimeSerie[i]
                Velocidad[AvalNumber,2] = nflag
                TiempoAval[AvalNumber,2] = nflag
                TiempoAval[AvalNumber,3] = 0
                XsecondAval = x1[i]
                YsecondAval = y1[i]
                TimeAval = TimeSerie[i]
            else
                if TimeSerie[i]== TimeAval
                    DistanciaAval = sqrt((x1[i]-XsecondAval)^2 + (y1[i]- YsecondAval)^2)
                    DeltaTiempo = .00000000000000000000000001
                    Velocidad[AvalNumber,1] = DistanciaAval/DeltaTiempo
                    TiempoAval[AvalNumber,1]= TimeSerie[i]
                    Velocidad[AvalNumber,2] = nflag
                    TiempoAval[AvalNumber,2] = nflag
                    TiempoAval[AvalNumber,3] = DeltaTiempo
                    XsecondAval = x1[i]
                    YsecondAval = y1[i]
                    TimeAval = TimeSerie[i] + .00000000000000000000000001       
                else
                    DistanciaAval = sqrt((x1[i]-XsecondAval)^2 + (y1[i]- YsecondAval)^2)
                    DeltaTiempo = TimeSerie[i] - TimeAval
                    Velocidad[AvalNumber,1] = DistanciaAval/DeltaTiempo
                    TiempoAval[AvalNumber,1] = TimeSerie[i]
                    Velocidad[AvalNumber,2] = nflag
                    TiempoAval[AvalNumber,2] = nflag
                    TiempoAval[AvalNumber,3] = DeltaTiempo
                    XsecondAval = x1[i]
                    YsecondAval = y1[i]
                    TimeAval = TimeSerie[i]
                end
            end
        end
    end

    contAva = 0
    contNorm = 0
    PyPlot.figure(7389,(10,10))
    for i=2:AvalNumber
        if TiempoAval[i,2] == 1
            cl = "or"
            nsi = 2
            if Velocidad[i,1]!=0
              contAva += 1
              PyPlot.plot(log10.(TiempoAval[i,1]), log10.(Velocidad[i,1]), cl, markersize=nsi,alpha=0.5,label = "Synthetic_quakes")
              PyPlot.hold("True")
              if contAva == 1
                  if ik == 1
                      PyPlot.legend(fontsize=14, markerscale=3)
                  end                 
              end
            end
        else    
            cl = "or"
            nsi = 2
            if Velocidad[i,1]!=0
              contNorm += 1
              PyPlot.plot(log10.(TiempoAval[i,1]), log10.(Velocidad[i,1]), cl, markersize=nsi,alpha=0.5)
              PyPlot.hold("True")
            end
        end        
    end
    PyPlot.plot(log10.(TiempoAval[:,1]), log10.(Velocidad[:,1]), "gray", linewidth=0.5,alpha=0.5)
    PyPlot.xlabel("log10(Time dimensionless)",fontsize = 17,weight="semibold")
    PyPlot.ylabel("log10(Inter-event rate)",fontsize = 17,weight="semibold")
    PyPlot.xticks(fontsize=15,weight="semibold")
    PyPlot.yticks(fontsize=15,weight="semibold")   
    PyPlot.savefig(pathresults*VecID*"-IntereventRate.pdf",dpi=200)
    PyPlot.hold("False")   
#   ====================================================================================================


#   Spatial distribution and Inter-event rate of simulated earthquakes 
    l=0
    NUM=0
    vecColor=zeros(smin,3)
    for i=1:smin
        vecColor[i,1]=rand()
        vecColor[i,2]=rand()
        vecColor[i,3]=rand()
    end
    ContFin=0
    VecSize= zeros(length(VecSpatialAva[:,1]))
    VecScatX = zeros(length(VecSpatialAva[:,1]))
    VecScatY = zeros(length(VecSpatialAva[:,1]))
    VecCoorEpicJ =zeros(length(VecSpatialAva[:,1]))
    VecCoorEpicI =zeros(length(VecSpatialAva[:,1]))
    VecTiem=zeros(length(VecSpatialAva[:,1]))
    VecAreaHB1 = zeros(length(VecSize))
    VecMagHB1 = zeros(length(VecSize))
    VecX=zeros(length(VecSize))
    VecY=zeros(length(VecSize))
    VecEpiI= zeros(length(VecSize))
    VecEpiJ= zeros(length(VecSize))
    VecSizeR=zeros(length(VecSize))
	  VecTime	=zeros(length(VecSize))
    VecVelocityAval	=zeros(length(VecSize),6)
    r=0

    if ik == 1
      PyPlot.figure(8839)
          for j=2:length(datos[:,1])
              if datos[j,5]==1
                  NUM=NUM+1
                  l=l+1
                  XX[l,3]=datos[j,3]
                  XX[l,1]=datos[j,6]
                  XX[l,2]=datos[j,7]
                  l1[l]=NUM
                  VecSpatialAva[NUM,1]=XX[l,3] #tiempo
                  VecSpatialAva[NUM,2]=XX[l,1] #x
                  VecSpatialAva[NUM,3]=XX[l,2] #y
              elseif datos[j,5]==0 && datos[j-1,5]==1
            # Plot las coordenadas espaciales de las AVALANCHAS en diferentes colores
                  ContFin+=1
                  VecScatX[ContFin]=XX[1,1]
                  VecScatY[ContFin]=XX[1,2]
                  VecTiem[ContFin]=XX[1,3]
                  VecSize[ContFin]=l
                  AreaHB1 = CellSize # Area of a unitary cell 
                  PyPlot.grid("on")
                  if VecSize[ContFin]>=1
                      r=r+1
                      VecAreaHB1[r] = AreaHB1 * VecSize[ContFin]
  #                   Ramirez relation (Rodriguez Q and Otemöller, 2013)
                      VecMagHB1[r] = (2/3 * log10( (VecAreaHB1[r]/(7.78*1.0e-9)).^(1/0.550)) ) - 6.07
  #                   
                      if 2.0<VecMagHB1[r]<=2.5
                          A="y"
                      elseif 2.5<VecMagHB1[r]<=3
                          A="b"
                      elseif 3<VecMagHB1[r]<=4
                          A="r"
                      elseif 4<VecMagHB1[r]<=5
                          A="g"
                      elseif 5<VecMagHB1[r]
                          A="m"
                      end
                  VecX[r] = VecScatX[ContFin]
                  VecY[r] = VecScatY[ContFin]
                  VecSizeR[r]= VecSize[ContFin]
                  VecTime[r]= VecTiem[ContFin]

  #                 PyPlot.subplot(211)
                  
                  if VecMagHB1[r] > 1
                      PyPlot.plot(XX[1:l,1],XX[1:l,2], linestyle="None", "o", markersize=3, color=vecColor[j,:],alpha=0.5)
                      PyPlot.xlabel("Nx [cells]",fontsize=16,weight="bold")
                      PyPlot.ylabel("Ny [cells]",fontsize=16,weight="bold")
                  end    
                  
  #                 PyPlot.subplot(212)
  #                 l1Int1 = round(Int,l1[1])
  #                 l1Intl = round(Int,l1[l])
  #                 PyPlot.plot(log10.(TiempoAval[l1Int1:l1Intl,1]),log10.(Velocidad[l1Int1:l1Intl,1]), "ok", markersize=1,alpha=0.5)
  #                 PyPlot.xlabel("log10(Time[Dimensionless])",fontsize=16,weight="bold")
  #                 PyPlot.ylabel("log10(Inter-event rate)",fontsize=16,weight="bold")   
  #                 PyPlot.hold("True")

              end

              XX[1:l,:]=0.0
              l=0
              l1[1:l]=0.0
          end
      end
      
        xticks(fontsize=14,weight="demi")
        yticks(fontsize=14,weight="demi")
        PyPlot.savefig(pathresults*VecID*"-SpatialDistribution.pdf",dpi=200)
        PyPlot.hold("False")
    end
 # ====================================================================================================       

end
