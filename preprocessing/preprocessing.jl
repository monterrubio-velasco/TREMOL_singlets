# Author:  Marisol Monterrubio-Velasco
# Contact: marisol.monterrubio@bsc.es
# GNU GENERAL PUBLIC LICENSE
# ==========================
# Version 3, 29 June 2007
# Copyright © 2007 Free Software Foundation, Inc. <http://fsf.org/>

# preprocessing.jl : Define the input values and its functions and pass it to the TREMOL_main.jl script

#Import the modules to plot the results and the configuration
#   using PyPlot
#   using PyCall
#   @pyimport matplotlib.colors as matcolors
#   @pyimport matplotlib as mpl
#   @pyimport matplotlib.patches as patches
# 
#Import function
#- include(joinpath("/YOURPATH/TREMOL_singlets_SUB3/main/","TREMOL_Singlets.jl"))
# 
#Parameters
# =======================================================================
# Parameters to assigned the initial conditions to the simulated domain 
# ========================================================================
#   - `VecID ::String: Name to identify the simulated event.
#   - `fhi_asp::Float`: percentage of load to be transferred for any ruptured cell in the asperity domain [0,1]
#   - `fhi_bkg::Float or Vector`: percentage of load to be transferred for any ruptured cell in the background [0,1]
#   - `strength_asp::Integer`: value to define the strength in the asperity cells >= 1.
#   - `nbox::Integer`: Number of cells assigned to each Seismic Source.
#   - `nk::Integer`: Number of times that will execute a same experiment to obtain statistic 
#   - `RecFact::Float`: Factor that converts asperity and background area from a square shape Ra=1 to a rectangular shape for Ra>1 
# ====================================================
# Parameters coming from the finite fault source method
# ====================================================
#   - `DurTeo::Float Vector`: rupture duration given by finite fault computation [seconds]
#   - `VelAsp :: Float Vector`: rupture velocity given by finite fault computation [km/s]
#   - `Leff::Float Vector`: Effective length size of the effective rupture area [km]
#   - `Weff::Float Vector`: Effective wide size of the effective rupture area [km]
#   - `VectorSeff::Float Vector`: Effective area [km²]
#   - `SaOri::Float Vector`: Ratio of the asperity size
#   - `YsizeOri::Float Vector`: Asperity area [km²]
# ====================================================

#...
# """

YOURPATH = ""   #Here define your path

using PyPlot
using PyCall
@pyimport matplotlib.colors as matcolors
@pyimport matplotlib as mpl
@pyimport matplotlib.patches as patches
include(joinpath("YOURPATH/TREMOL_singlets_SUB3/main/","TREMOL_Singlets.jl"))
include(joinpath("YOURPATH/TREMOL_singlets_SUB3/resultsPostProcessing/","eventsRealSyntheticHistogramStatistics.jl"))

# Input values 
# VecID = String Array # ["An ID name"]  
# fhi_bkg = Float # 0.67
# fhi_asp = Float # 0.90
# strength_asp = Integer # 4
# nbox = Integer # 100
# nk = Integer # 1
# DurTeo = Float Array # [40] 
# VelAsp = Float Array # [3.2]
# Leff = Float Array # [34.47] 
# Weff = Float Array # [17.81]   
# VectorSeff = Float Array # [Leff[1].*Weff[1]]   
# SaOri = Float Array # [0.23]
# YsizeOri = Float Array #[Leff[1].*Weff[1].*SaOri[1]]  
# RecFact = Float

# ==========================================================
# # An example ...
# =========================================================

VecID = "28" #  --> 20/03/2012 M = 7.4, 
fhi_bkg =  0.67
fhi_asp =  0.90
strength_asp = 4
nbox = 100
nk = 1
DurTeo =  [30] 
VelAsp = [2.7]
Leff = [54.94] 
Weff = [53.59]   
VectorSeff = [Leff[1].*Weff[1]]   
SaOri = [0.26]
YsizeOri = [Leff[1].*Weff[1].*SaOri[1]]  
RecFact = 1.0

MatrizMagnitudesHistoSomer,MatrizMagnitudesHistoMai,MatrizMagnitudesHistoMai_VL,MatrizMagnitudesHistoRami = TREMOL_Singlets(VecID,fhi_asp,fhi_bkg,strength_asp,DurTeo,VelAsp,Leff,Weff,VectorSeff,SaOri,YsizeOri,nbox,nk,RecFact)
  
# Writing the synthetic results in the resultsPostProcessinn folder for further analysis

writedlm(joinpath("YOURPATH/TREMOL_singlets_SUB3/resultsPostProcessing/", "Cels"*string(nbox)*"_MagniHistoSomer-Aeff"*"-"*VecID*"-SUB3.dat"), MatrizMagnitudesHistoSomer)
writedlm(joinpath("YOURPATH/TREMOL_singlets_SUB3/resultsPostProcessing/", "Cels"*string(nbox)*"_MagniHistoMai_L-Aeff"*"-"*VecID*"-SUB3.dat"), MatrizMagnitudesHistoMai)
writedlm(joinpath("YOURPATH/TREMOL_singlets_SUB3/resultsPostProcessing/", "Cels"*string(nbox)*"_MagniHistoMai_VL-Aeff"*"-"*VecID*"-SUB3.dat"), MatrizMagnitudesHistoMai_VL)
writedlm(joinpath("YOURPATH/TREMOL_singlets_SUB3/resultsPostProcessing/", "Cels"*string(nbox)*"_MagniHistoRami-Aeff"*"-"*VecID*"-SUB3.dat"), MatrizMagnitudesHistoRami)
