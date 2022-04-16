from optparse import OptionParser

# plotting :
from pylab import *
import numpy
import math
import matplotlib.pyplot as plt
import matplotlib.dates as md  # to convert unix time to date on X-axis 
import os
import copy
import time
import math
from scipy.interpolate import interp1d # for interpolation 

import matplotlib
matplotlib.rc('xtick', labelsize=25) 
matplotlib.rc('ytick', labelsize=25) 

# import radec2azim

# time operations :
from astropy.coordinates import SkyCoord, EarthLocation
from astropy.time import Time
from datetime import datetime
MWA_POS=EarthLocation.from_geodetic(lon="116:40:14.93",lat="-26:42:11.95",height=377.8)

 
def parse_options(idx):
   global debug_level

   usage="Usage: %prog [options]\n"
   usage += "Plotting FRB rates using results of different papers scaled down to 150 or 200 MHz :\n"
   
   parser = OptionParser(usage=usage,version=1.00)

   # General parameters :
#   parser.add_option('-s','--station_name','--station',dest="station_name",default="AAVS2", help="Station name [default %default]")

   # TEST : azim=0, za=0, lst=15.4 
#   parser.add_option('-p','--plot','--do_plot','--do_plots',action="store_true",dest="do_plot",default=False, help="Plot")
#   parser.add_option('-a','--azim','--azim_deg',dest="azim_deg",default=None, help="Pointing direction azimuth in degrees [default %default]",metavar="float",type="float")
   parser.add_option('-f','--freq','--freq_mhz',dest="freq_mhz",default=154.00, help="Survey frequency in MHz [default %default MHz]",metavar="float",type="float")
   (options, args) = parser.parse_args(sys.argv[idx:])

   print("#####################################################################################################")
   print("PARAMETERS:")
   print("#####################################################################################################")
   print("Frequency = %.2f [MHz]" % (options.freq_mhz))
   print("#####################################################################################################")
   
   return (options, args)



#####################################################################################################

#####################################################################################################
def frb_rate_euclidean( fluence, fluence_ref, frb_rate_ref , debug=True ) :
   if debug : 
      print("DEBUG : frb_rate_euclidean(%.4f,%.4f,%.4f)" % (fluence, fluence_ref, frb_rate_ref))
#   frb_rate = frb_rate_ref*math.pow( fluence/fluence_ref , -1.5 )
#   frb_rate = 1
   frb_rate = float(frb_rate_ref)*math.pow( (float(fluence)/float(fluence_ref)) , -1.5 )
   
   if debug : 
      print("DEBUG :    %.8f " % (frb_rate))
   
   return frb_rate

#####################################################################################################
# Given certain FRB rate at, for example, 700 FRBs / day / sky above fluence 700 Jy ms , it returns 
# FRB rates for various fluences in the range min_fluence to max_fluence 
#####################################################################################################
def calc_frb_rates_vs_fluence( freq_target,                    # target frequency
                               fluence_ref, frb_rate_ref ,     # reference fluence and FRB rate 
                               freq_ref=154, spectral_index=0, # reference frequency and spectral index for scaling, 0 - flat spectrum 
                               min_fluence=1, max_fluence=10000, fluence_step=1 # range to be outputted 
                             ): 
   out_fluence = []
   out_rates   = []

   fluence_freq = fluence_ref*math.pow( (float(freq_target)/float(freq_ref)) , spectral_index )
   
   fluence = min_fluence
   while fluence < max_fluence :
      frb_rate = frb_rate_euclidean( fluence, fluence_freq, frb_rate_ref )
      
      out_fluence.append( fluence )
      out_rates.append( frb_rate )
   
      fluence += fluence_step   

      if fluence > 10 :
         fluence_step = 10

      if fluence > 100 :
         fluence_step = 100

      if fluence > 1000 :
         fluence_step = 1000
         

   return (out_fluence,out_rates)   
   

def do_plots(options) :
   fig_size_x=20
   fig_size_y=10
   legend_list = []
   legend_location = "upper right"

   plt.figure( figsize=( fig_size_x , fig_size_y ) )
   fig = plt.gcf()
   
   # add points from Tingay et al (2015) :
   fluence_limit_tingay = 700
   frb_rate_tingay = 700
   point_tingay = 'v'
   (fluence_tingay,frb_rates_tingay) = calc_frb_rates_vs_fluence( options.freq_mhz, fluence_limit_tingay, frb_rate_tingay, spectral_index=0 )
   ax =  plt.plot( fluence_tingay,frb_rates_tingay, linestyle='dashed', color='green', linewidth=1 )   
   ax =  plt.plot( fluence_tingay,frb_rates_tingay, point_tingay, color='green', linewidth=2, markersize=5 )
   legend_list.append("Tingay et al.,2015, 154 MHz")
   legend_list.append("Tingay et al.,2015, 154 MHz")
   plt.plot( [fluence_limit_tingay] , [frb_rate_tingay] , point_tingay , color='green', markersize=12 )
   legend_list.append("Tingay et al.,2015, 154 MHz")
   
   # Shannon et al. 2018 
   # 37 +/- 8 /day/sky 
   # above threshold = 26 Jy ms (w/1.26 ms)^{-1/2}
   # Frequency spectral index scaling_index=-2.1
   fluence_limit_rshannon = 26
   frb_rate_rshannon = 37
   spectral_index_rshannon = -2.1
   freq_shannon = 1400.00
   point_rshannon = '+'

   # using Ryan's spectral index from ASKAP :
   (fluence_rshannon,frb_rates_rshannon) = calc_frb_rates_vs_fluence( options.freq_mhz, fluence_limit_rshannon, frb_rate_rshannon, spectral_index=spectral_index_rshannon, freq_ref=freq_shannon )
   ax_rshannon =  plt.plot( fluence_rshannon,frb_rates_rshannon, linestyle='dashed', color='orange', linewidth=2, markersize=12 )
   legend_list.append("Shannon et al.,2018, #alpha=-2.1")
#   plt.plot( [fluence_limit_rshannon] , [frb_rate_rshannon] , point_rshannon , color='orange', markersize=20 )

   # using spectral index = -1
   (fluence_rshannonM1,frb_rates_rshannonM1) = calc_frb_rates_vs_fluence( options.freq_mhz, fluence_limit_rshannon, frb_rate_rshannon, spectral_index=-1 , freq_ref=freq_shannon )
   ax_rshannonM1 =  plt.plot( fluence_rshannonM1,frb_rates_rshannonM1, linestyle='-.', color='orange', linewidth=2, markersize=12 )
   legend_list.append("Shannon et al.,2018, #alpha=-1")

   # using flat (0) spectral index :
   (fluence_rshannon0,frb_rates_rshannon0) = calc_frb_rates_vs_fluence( options.freq_mhz, fluence_limit_rshannon, frb_rate_rshannon, spectral_index=0 , freq_ref=freq_shannon )
   ax_rshannon0 =  plt.plot( fluence_rshannon0,frb_rates_rshannon0, linestyle='dotted', color='orange', linewidth=2, markersize=12 )
   legend_list.append("Shannon et al.,2018, #alpha=0")


   plt.legend( legend_list, loc=legend_location, fontsize=20)

   
   
   plt.xlabel('Fluence [Jy ms]' , fontsize=30 )
   plt.ylabel('#FRBs / day / sky' , fontsize=30 )
   plt.xscale('log')
   plt.yscale('log')

   plt.savefig( "frb_rates.png" )
   plt.show()   



if __name__ == "__main__":
    (options, args) = parse_options(1)
    
    do_plots(options)
    
        
