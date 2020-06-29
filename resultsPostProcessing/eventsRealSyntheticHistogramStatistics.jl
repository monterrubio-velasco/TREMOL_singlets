#Extraction of the mean and standard deviation of the histograms of magnitude
using PyPlot

function eventsRealSyntheticHistogramStatistics(YOURPATH,EventID)

    Y0 = zeros(23,65);
    Y1 = zeros(23,65);
    Y2 = zeros(23,65);
    Y3 = zeros(23,65);


    Event = "Evento22"

    X0 = readdlm("YOURPATH/TREMOL_singlets_SUB3/resultsPostProcessing/Synthetic_DataBase/Cels200_MagniHistoMai_L_Aeff-MesSa-"*Event*"-SUB3-2.dat")
    X1 = readdlm("YOURPATH/TREMOL_singlets_SUB3/resultsPostProcessing/Synthetic_DataBase/Cels200_MagniHistoMai_VL_Aeff-MesSa-"*Event*"-SUB3-2.dat")
    X2 = readdlm("YOURPATH/TREMOL_singlets_SUB3/resultsPostProcessing/Synthetic_DataBase/Cels200_MagniHistoSomerAeff-MesSa-"*Event*"-SUB3-2.dat")
    X3 = readdlm("YOURPATH/TREMOL_singlets_SUB3/resultsPostProcessing/Synthetic_DataBase/Cels200_MagniHistoRamiAeff-MesSa-"*Event*"-SUB3-2.dat")

    for i = 1:65
      Y0[1,i] = mean(X0[:,i])
      Y0[2,i] = std(X0[:,i])
      Y1[1,i] = mean(X1[:,i])
      Y1[2,i] = std(X1[:,i])
      Y2[1,i] = mean(X2[:,i])
      Y2[2,i] = std(X2[:,i])
      Y3[1,i] = mean(X3[:,i])
      Y3[2,i] = std(X3[:,i])  
    end

    Event = "Evento28"

    X0 = readdlm("YOURPATH/TREMOL_singlets_SUB3/resultsPostProcessing/Synthetic_DataBase/Cels200_MagniHistoMai_L_Aeff-MesSa-"*Event*"-SUB3-1.dat")
    X1 = readdlm("YOURPATH/TREMOL_singlets_SUB3/resultsPostProcessing/Synthetic_DataBase/Cels200_MagniHistoMai_VL_Aeff-MesSa-"*Event*"-SUB3-1.dat")
    X2 = readdlm("YOURPATH/TREMOL_singlets_SUB3/resultsPostProcessing/Synthetic_DataBase/Cels200_MagniHistoSomerAeff-MesSa-"*Event*"-SUB3-1.dat")
    X3 = readdlm("YOURPATH/TREMOL_singlets_SUB3/resultsPostProcessing/Synthetic_DataBase/Cels200_MagniHistoRamiAeff-MesSa-"*Event*"-SUB3-1.dat")

    Y0[3,:] = X0[21,:]
    Y0[4,:] = X0[22,:]
    Y1[3,:] = X1[21,:]
    Y1[4,:] = X1[22,:]
    Y2[3,:] = X2[21,:]
    Y2[4,:] = X2[22,:]
    Y3[3,:] = X3[21,:]
    Y3[4,:] = X3[22,:]  

    Event = "Evento24"

    X0 = readdlm("YOURPATH/TREMOL_singlets_SUB3/resultsPostProcessing/Synthetic_DataBase/Cels200_MagniHistoMai_L_Aeff-MesSa-"*Event*"-SUB3-1.dat")
    X1 = readdlm("YOURPATH/TREMOL_singlets_SUB3/resultsPostProcessing/Synthetic_DataBase/Cels200_MagniHistoMai_VL_Aeff-MesSa-"*Event*"-SUB3-1.dat")
    X2 = readdlm("YOURPATH/TREMOL_singlets_SUB3/resultsPostProcessing/Synthetic_DataBase/Cels200_MagniHistoSomerAeff-MesSa-"*Event*"-SUB3-1.dat")
    X3 = readdlm("YOURPATH/TREMOL_singlets_SUB3/resultsPostProcessing/Synthetic_DataBase/Cels200_MagniHistoRamiAeff-MesSa-"*Event*"-SUB3-1.dat")

    Y0[5,:] = X0[21,:]
    Y0[6,:] = X0[22,:]
    Y1[5,:] = X1[21,:]
    Y1[6,:] = X1[22,:]
    Y2[5,:] = X2[21,:]
    Y2[6,:] = X2[22,:]
    Y3[5,:] = X3[21,:]
    Y3[6,:] = X3[22,:]
    
    
    Y0[9,:] = X0[23,:]
    Y1[9,:] = X1[23,:]
    Y2[9,:] = X2[23,:]
    Y3[9,:] = X3[23,:]

    Event = "Evento20"

    X0 = readdlm("YOURPATH/TREMOL_singlets_SUB3/resultsPostProcessing/Synthetic_DataBase/Cels200_MagniHistoMai_L_Aeff-MesSa-"*Event*"-SUB3-2.dat")
    X1 = readdlm("YOURPATH/TREMOL_singlets_SUB3/resultsPostProcessing/Synthetic_DataBase/Cels200_MagniHistoMai_VL_Aeff-MesSa-"*Event*"-SUB3-2.dat")
    X2 = readdlm("YOURPATH/TREMOL_singlets_SUB3/resultsPostProcessing/Synthetic_DataBase/Cels200_MagniHistoSomerAeff-MesSa-"*Event*"-SUB3-2.dat")
    X3 = readdlm("YOURPATH/TREMOL_singlets_SUB3/resultsPostProcessing/Synthetic_DataBase/Cels200_MagniHistoRamiAeff-MesSa-"*Event*"-SUB3-2.dat")

    for i = 1:65
      Y0[7,i] = mean(X0[:,i])
      Y0[8,i] = std(X0[:,i])
      Y1[7,i] = mean(X1[:,i])
      Y1[8,i] = std(X1[:,i])
      Y2[7,i] = mean(X2[:,i])
      Y2[8,i] = std(X2[:,i])
      Y3[7,i] = mean(X3[:,i])
      Y3[8,i] = std(X3[:,i])  
    end
    
    RealData = zeros(10)
    
    for i = 1:65
    
      if EventID == "20"
      
        Y0[11,i] = Y0[7,i]  #suma de los valores medios Y0[1,i] #
        Y0[12,i] = Y0[8,i]
        Y1[11,i] = Y1[7,i]
        Y1[12,i] = Y1[8,i]
        Y2[11,i] = Y2[7,i]
        Y2[12,i] = Y2[8,i]
        Y3[11,i] = Y3[7,i]
        Y3[12,i] = Y3[8,i]
        RealData = readdlm("YOURPATH/TREMOL_singlets_SUB3/resultsPostProcessing/Real_DataBase/Foreshocks-Event-20-MagniGR.dat")
        
      elseif EventID == "22"
      
        Y0[11,i] = Y0[1,i]  #suma de los valores medios Y0[1,i] #
        Y0[12,i] = Y0[2,i]
        Y1[11,i] = Y1[1,i]
        Y1[12,i] = Y1[2,i]
        Y2[11,i] = Y2[1,i]
        Y2[12,i] = Y2[2,i]
        Y3[11,i] = Y3[1,i]
        Y3[12,i] = Y3[2,i]
        RealData = readdlm("YOURPATH/TREMOL_singlets_SUB3/resultsPostProcessing/Real_DataBase/Foreshocks-Depth_pm_8-Ev22--MagniGR.dat")
      
      elseif EventID == "24"
      
        Y0[11,i] = Y0[5,i]  #suma de los valores medios Y0[1,i] #
        Y0[12,i] = Y0[6,i]
        Y1[11,i] = Y1[5,i]
        Y1[12,i] = Y1[6,i]
        Y2[11,i] = Y2[5,i]
        Y2[12,i] = Y2[6,i]
        Y3[11,i] = Y3[5,i]
        Y3[12,i] = Y3[6,i]
        RealData = readdlm("YOURPATH/TREMOL_singlets_SUB3/resultsPostProcessing/Real_DataBase/Foreshocks-Depth_pm_8-Ev24--MagniGR.dat")
        
      elseif EventID == "28"
      
        Y0[11,i] = Y0[3,i]  #suma de los valores medios Y0[1,i] #
        Y0[12,i] = Y0[4,i]
        Y1[11,i] = Y1[3,i]
        Y1[12,i] = Y1[4,i]
        Y2[11,i] = Y2[3,i]
        Y2[12,i] = Y2[4,i]
        Y3[11,i] = Y3[3,i]
        Y3[12,i] = Y3[4,i]
        RealData = readdlm("YOURPATH/TREMOL_singlets_SUB3/resultsPostProcessing/Real_DataBase/Foreshocks-Depth_pm_8-Ev28--MagniGR.dat")
        
      
      elseif EventID == "All"
      
        Y0[11,i] = Y0[1,i]+Y0[3,i]+Y0[5,i]+Y0[7,i] #suma de los valores medios Y0[1,i] #
        Y0[12,i] = Y0[2,i]+Y0[4,i]+Y0[6,i]+Y0[8,i]#suma de errores Y0[2,i]  #
        Y1[11,i] = Y1[1,i]+Y1[3,i]+Y1[5,i]+Y1[7,i]
        Y1[12,i] = Y1[2,i]+Y1[4,i]+Y1[6,i]+Y1[8,i]
        Y2[11,i] = Y2[1,i]+Y2[3,i]+Y2[5,i]+Y2[7,i]
        Y2[12,i] = Y2[2,i]+Y2[4,i]+Y2[6,i]+Y2[8,i]
        Y3[11,i] = Y3[1,i]+Y3[3,i]+Y3[5,i]+Y3[7,i]
        Y3[12,i] = Y3[2,i]+Y3[4,i]+Y3[6,i]+Y3[8,i]
        RealData = readdlm("YOURPATH/TREMOL_singlets_SUB3/resultsPostProcessing/Real_DataBase/4Foreshocks-Depth_pm_8-MagniGR.dat")
      end
      
    end


    sum0 = 0
    sumerr0 =0 
    sum1=0
    sumerr1=0
    sum2=0
    sumerr2=0
    sum3=0
    sumerr3=0

    for i = 1:65
            
        sum0 = sum(Y0[11,i:65])
        Y0[13,i] = sum0
        sumerr0 = sum(Y0[12,i:65])
        Y0[14,i] = sumerr0
        
        sum1 = sum(Y1[11,i:65])
        Y1[13,i] = sum1
        sumerr1 = sum(Y1[12,i:65])
        Y1[14,i] = sumerr1    
        
        sum2 = sum(Y2[11,i:65])
        Y2[13,i] = sum2
        sumerr2 = sum(Y2[12,i:65])
        Y2[14,i] = sumerr2
            
        sum3 = sum(Y3[11,i:65])
        Y3[13,i] = sum3
        sumerr3 = sum(Y3[12,i:65])
        Y3[14,i] = sumerr3
    
      
    end


    Y0[15,:] = Y0[13,:] + (Y0[14,:]/2)
    Y0[16,:] = Y0[13,:] - (Y0[14,:]/2)
    Y1[15,:] = Y1[13,:] + (Y1[14,:]/2)
    Y1[16,:] = Y1[13,:] - (Y1[14,:]/2)
    Y2[15,:] = Y2[13,:] + (Y2[14,:]/2)
    Y2[16,:] = Y2[13,:] - (Y2[14,:]/2)
    Y3[15,:] = Y3[13,:] + (Y3[14,:]/2)
    Y3[16,:] = Y3[13,:] - (Y3[14,:]/2)

    

    #  Event = "28"
    # #############################################
    # #NonAsperity_Data
    # ###############################################
    # Y0non = readdlm("Event28HistoMagni_Mai_L-NonAsp.dat")
    # Y1non = readdlm("Event28HistoMagni_Mai_VL-NonAsp.dat")
    # Y2non = readdlm("Event28HistoMagni_Somer-NonAsp.dat")
    # Y3non = readdlm("Event28HistoMagni_Rami-NonAsp.dat")
    # ###############################################
    # 
    # #NonStrength_Data
    # ###############################################
    # Y0nonStr = readdlm("EventNONStrength-28-HistoMagni_Mai_L.dat")
    # Y1nonStr = readdlm("EventNONStrength-28-HistoMagni_Mai_VL.dat")
    # Y2nonStr = readdlm("EventNONStrength-28-HistoMagni_Somer.dat")
    # Y3nonStr = readdlm("EventNONStrength-28-HistoMagni_Rami.dat")
    ###############################################


    ########################################################################################################################3
    ###   Mai et al, 2005
    #####################################################################
    ini = 2
    fin = 65
    PyPlot.figure(738,(7,6))
    # ylim(0.5,1000)
    xlim(2,9)
    PyPlot.hold("on")

    PyPlot.semilogy(Y0[9,:],RealData[1:65],label="real data",linewidth=2.0,"ob--",ms=5,alpha=0.5)
    PyPlot.semilogy(Y0[9,ini:fin],Y0[13,ini:fin],label="synt. data (Eq. 3)",linewidth=2.5,"sk-",ms=3) 
    PyPlot.semilogy(Y0[9,ini:fin],Y0[15,ini:fin],label=L"$mean \pm \sigma$",linewidth=2.0,"-.k",ms=3)
    PyPlot.semilogy(Y0[9,ini:fin],Y0[16,ini:fin],linewidth=2.0,"-.k") #,label=L"$mean - \sigma$",linewidth=1.0,":>k",ms=3)
   
    xticks(fontsize=14,weight="semibold")
    yticks(fontsize=14,weight="semibold")
    PyPlot.xlabel("Magnitude",fontsize=14,weight="semibold")
    PyPlot.ylabel("Cumulative Number",fontsize=14,weight="semibold")
    legend(fontsize=14)

     PyPlot.savefig("YOURPATH/TREMOL_singlets_SUB3/resultsPostProcessing/Foreshocks-StatdifSa-Event"*EventID*"-SUB3-N200_Mai-L.pdf",dpi = 400)

    ###################################################################################################################################### Mai et al, 2005
    ######################################################################################################################################

  
    PyPlot.figure(740,(7,6))
    # ylim(Y1[19,52],6000)
#     ylim(0.5,1000)
    xlim(2,9)
    PyPlot.hold("on")

    PyPlot.semilogy(Y0[9,:],RealData[1:65],label="real data",linewidth=2.0,"ob--",ms=5,alpha=0.5)
    PyPlot.semilogy(Y1[9,ini:fin],Y1[13,ini:fin],label="synt. data (Eq. 4)",linewidth=2.5,"sk-",ms=3)
    PyPlot.semilogy(Y1[9,ini:fin],Y1[15,ini:fin],label=L"$mean \pm \sigma$",linewidth=2.0,"-.k",ms=3)
    PyPlot.semilogy(Y1[9,ini:fin],Y1[16,ini:fin],linewidth=2.0,"-.k") #,label=L"$mean - \sigma$",linewidth=1.0,":>k",ms=3)
   

    xticks(fontsize=14,weight="semibold")
    yticks(fontsize=14,weight="semibold")
    PyPlot.xlabel("Magnitude",fontsize=14,weight="semibold")
    PyPlot.ylabel("Cumulative Number",fontsize=14,weight="semibold")
    legend(fontsize=14)
    PyPlot.savefig("YOURPATH/TREMOL_singlets_SUB3/resultsPostProcessing/Foreshocks-StatdifSa-Event"*EventID*"-SUB3-N200_Mai-VL.pdf",dpi = 400)

    ######################################################################################################################################  Somerville et al., 1999a
    ######################################################################################################################################


    # PyPlot.figure(741,(20,10))
    PyPlot.figure(741,(7,6))
    # ylim(Y2[19,finplot],6000)
#     ylim(0.5,1000)
    xlim(2,9)
    PyPlot.hold("on")

    PyPlot.semilogy(Y0[9,:],RealData[1:65],label="real data",linewidth=2.0,"ob--",ms=5,alpha=0.5)
    PyPlot.semilogy(Y2[9,ini:fin],Y2[13,ini:fin],label="synt. data (Eq. 2)",linewidth=2.5,"sk-",ms=3)
    PyPlot.semilogy(Y2[9,ini:fin],Y2[15,ini:fin],label=L"$mean \pm \sigma$",linewidth=2.0,"-.k",ms=3)
    PyPlot.semilogy(Y2[9,ini:fin],Y2[16,ini:fin],linewidth=2.0,"-.k") #,label=L"$mean - \sigma$",linewidth=1.0,":>k",ms=3)
   
    xticks(fontsize=14,weight="semibold")
    yticks(fontsize=14,weight="semibold")
    PyPlot.xlabel("Magnitude",fontsize=14,weight="semibold")
    PyPlot.ylabel("Cumulative Number",fontsize=14,weight="semibold")
    legend(fontsize=14)
    PyPlot.savefig("YOURPATH/TREMOL_singlets_SUB3/resultsPostProcessing/Foreshocks-StatdifSa-Event"*EventID*"-SUB3-N200_Somer.pdf",dpi = 400)

    ###################################################################################################################################           Ramirez-Gaytan et al. 2014, magnitude area relation
    ######################################################################################################################################

    PyPlot.figure(742,(7,6))
#     ylim(0.01,10000)
    xlim(2,9)

    PyPlot.semilogy(Y0[9,:],RealData[1:65],label="real data",linewidth=2.0,"ob--",ms=5,alpha=0.5)
    PyPlot.semilogy(Y3[9,ini:fin],Y3[13,ini:fin],linewidth=2.5,"sk-",ms=3, label = "synth. data Asp.")
    PyPlot.semilogy(Y3[9,ini:fin],Y3[15,ini:fin],label=L"$mean \pm \sigma$",linewidth=2.0,"-.k",ms=3)
    PyPlot.semilogy(Y3[9,ini:fin],Y3[16,ini:fin],linewidth=2.0,"-.k")

    xticks(fontsize=14,weight="semibold")
    yticks(fontsize=14,weight="semibold")
    PyPlot.xlabel("Magnitude",fontsize=14,weight="semibold")
    PyPlot.ylabel("Cumulative Number",fontsize=14,weight="semibold")
    legend(fontsize=14)
    PyPlot.savefig(joinpath("YOURPATH/TREMOL_singlets_SUB3/resultsPostProcessing/", "Foreshocks-Event"*EventID*"-SUB3-N200_Ramirez.pdf"))
    # 
    # ###############################################################################################################
    # ##############################################################################################################
end
