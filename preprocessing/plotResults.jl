# Author:  Marisol Monterrubio-Velasco
# Contact: marisol.monterrubio@bsc.es
# GNU GENERAL PUBLIC LICENSE
# ==========================
# Version 3, 29 June 2007
# Copyright © 2007 Free Software Foundation, Inc. <http://fsf.org/>

# plotResults.jl : Plot the results obtained in TREMOL for the SUB3 analysis 

#Import the modules to plot the results and the configuration
#   using PyPlot
#   using PyCall
#   @pyimport matplotlib.colors as matcolors
#   @pyimport matplotlib as mpl
#   @pyimport matplotlib.patches as patches
# 
#Import function
#- include(joinpath("/YOURPATH/TREMOL_singlets_SUB3/resultsPostProcessing/","eventsRealSyntheticHistogramStatistics.jl"))
# 
#Parameters
# =======================================================================
# Parameters to assigned the initial conditions to the simulated domain 
# ========================================================================
# `EventID  :: String: Name that defines the event that is analyzed
# `YOURPATH :: String: Path to results folder and script
#
# ====================================================
# Return histogram frequency results for differen magnitude-area relations
# ====================================================
#...
# """

YOURPATH = " "  #Here define your path

using PyPlot
using PyCall
@pyimport matplotlib.colors as matcolors
@pyimport matplotlib as mpl
@pyimport matplotlib.patches as patches
include(joinpath("/home/mmonterr/Documentos/TREMOL/Tremol_SUB3Paper-GMD/code/TREMOL_singlets_SUB3/resultsPostProcessing/","eventsRealSyntheticHistogramStatistics.jl"))

# Choose the event to plot the ID:
# EventID = "20"  # --> 14/09/1995 M = 7.4,
# EventID = "22"  # --> 25/02/1996 M = 7.1,
# EventID = "24"  #--> 19/07/1997 M = 6.5,
# EventID = "28"  #--> 20/03/2012 M = 7.4,
EventID = "All" # --> cumulative sum of 4 events

eventsRealSyntheticHistogramStatistics(YOURPATH,EventID)

    
