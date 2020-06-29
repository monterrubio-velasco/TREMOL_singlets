# Author:  Marisol Monterrubio-Velasco
# Contact: marisol.monterrubio@bsc.es
# GNU GENERAL PUBLIC LICENSE
# ==========================
# Version 3, 29 June 2007
# Copyright © 2007 Free Software Foundation, Inc. <http://fsf.org/>

function gutenberRichAspDef(ik,VecMagHB1,nflag,AreaSUB,pathresults)

# """
# gutenberRichAspDef.jl: Function that computes the Gutenberg-Richter fit relation using the method of least squares 
 
# Parameters
# - `VecMagHB1::Array`: Vector that contains the equivalent magnitudes
# - `nflag::Integer`: Identifier
# - `AreaSUB::Float`: Cell area ion km²

# return
# - `paramGRAreaWY::Array`: results of the Gutenber-Richter fitting using the least square method
  ##1. paramGRAreaWY[1]=magMax
  ##2. paramGRAreaWY[2]=magMin
  ##3. paramGRAreaWY[3]=p1 (a-value)
  ##4. paramGRAreaWY[4]=p2 (b-value)
  ##5. paramGRAreaWY[5]=rho (correlation coeficient)
# - `cuenrep::Integer`: number of times that a same frequency of magnitudes is repeated for different magnitudes  
# - `vecWrite::Array`: vector that contains the frequency-magnitude data.
# """

    VecMagni=copy(VecMagHB1)
    contador=0
    paramGRAreaWY=zeros(5)
    srand(1)
    colorVect=rand(25, 3)
    magMax=maximum(VecMagni)
    magMin=minimum(VecMagni)
    magRang=magMax-magMin
    numMag=round(Int,magRang/0.25)
    limMin=zeros(numMag)
    NumSismos=zeros(numMag)
    vecWrite=zeros(numMag,2)

    for i=1:numMag
        limMin[i]=magMin+((i*0.25)-0.25)
        for j=1:length(VecMagni)
            if VecMagni[j]>=limMin[i]
                contador=contador+1
            end
        end
        NumSismos[i]=contador
        contador=0
    end
    cuenrep=0
    for i=2:numMag
        if NumSismos[i]==NumSismos[i-1]
            cuenrep += 1
        end
    end
    if magRang>0
        vecWrite[1:numMag,1]=limMin[1:numMag]
        vecWrite[1:numMag,2]=NumSismos[1:numMag]
        p1,p2 = polyfit(limMin[1:numMag],log10.(NumSismos[1:numMag]),1)
        rho=cor(limMin[1:numMag],log10.(NumSismos[1:numMag]))    
        
        if nflag == "Somer"
            cl = "g"
        elseif nflag == "Mai-L"
            cl = "b"
        elseif nflag == "Mai-VL"
            cl = "m"
        elseif nflag == "Ramirez"
            cl = "r"        
        end
       
        PyPlot.figure(77,(10,10))
        PyPlot.hold("on")
        PyPlot.plot(limMin[1:numMag],log10.(NumSismos[1:numMag]),label=nflag,linewidth=2.0,"o-",color = cl)
        if ik==1
            PyPlot.legend(fontsize=14, markerscale=2)
        end
        PyPlot.xlabel("Magnitude", fontsize = 17,weight="semibold")
        PyPlot.ylabel("log10(Cumulative number of events)",fontsize = 17,weight="semibold")
        PyPlot.xticks(fontsize=15,weight="semibold")
        PyPlot.yticks(fontsize=15,weight="semibold") 
        PyPlot.savefig(pathresults*"-GRfit.pdf",dpi = 400)
        
        paramGRAreaWY[1]=magMax
        paramGRAreaWY[2]=magMin
        paramGRAreaWY[3]=p1
        paramGRAreaWY[4]=p2
        paramGRAreaWY[5]=rho
    elseif magRang==0
        paramGRAreaWY[1]=0
        paramGRAreaWY[2]=0
        paramGRAreaWY[3]=0
        paramGRAreaWY[4]=0
        paramGRAreaWY[5]=0
    end

    return  paramGRAreaWY, cuenrep, vecWrite[1:numMag,:]
end


function polyfit(x, y, n)
    A = [ float(x[i])^p for i = 1:length(x), p = 0:n ]
    A \ y
end
