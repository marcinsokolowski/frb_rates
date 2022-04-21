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
  
# https://matplotlib.org/3.5.0/api/_as_gen/matplotlib.pyplot.legend.html   
# legend()
# legend(handles, labels)
# legend(handles=handles)
# legend(labels)
def do_plots(options) :
   fig_size_x=20
   fig_size_y=10
   legend_list = []
   plot_list   = []
   legend_location = "upper right"

   plt.figure( figsize=( fig_size_x , fig_size_y ) )
   fig = plt.gcf()
   
   ##################################################################################################################################################################       
   #                                   Tingay et al (2015) : /home/msok//Desktop/GRANTS/2020/LIEF/SECTIONS/references/LOW_FREQ_NON-DETECTIONS/Tingay_2015_AJ_150_199_MS.pdf   
   # Antonina R. : RFRB < 82 /sky/day at 182 MHz, above a fluence of F > 7980 Jy ms (not a strong limit -> not plotted)
   fluence_limit_tingay = 700
   frb_rate_tingay = 700
   point_tingay = 'v'
   (fluence_tingay,frb_rates_tingay) = calc_frb_rates_vs_fluence( options.freq_mhz, fluence_limit_tingay, frb_rate_tingay, spectral_index=0 )
   ax1 =  plt.plot( fluence_tingay,frb_rates_tingay, linestyle='dashed', color='green', linewidth=1 )   
#   ax2 =  plt.plot( fluence_tingay,frb_rates_tingay, point_tingay, color='green', linewidth=2, markersize=5 )
#   legend_list.append("Tingay et al.,2015, 154 MHz")
   ax3 = plt.plot( [fluence_limit_tingay] , [frb_rate_tingay] , point_tingay , linestyle='dashed', color='green', linewidth=2, markersize=10 )
   plot_list.append(ax3[0])   
   legend_list.append("Tingay et al.,2015, 154 MHz")
   plt.xlim(( 1, 10000 ))
#   ax1.set_xlim(1,10000)
#   legend_list.append("Tingay et al.,2015, 154 MHz")
   ##################################################################################################################################################################
   
   ##################################################################################################################################################################   
   #                             Shannon et al. 2018 page 5 in /home/msok/Desktop/GRANTS/2021/LIEF/SECTIONS/references/FRB/FRB_rates/ASKAP/Shannon_et_al_2018.pdf
   #
   # 37 +/- 8 /day/sky 
   # above threshold = 26 Jy ms (w/1.26 ms)^{-1/2}
   # Frequency spectral index scaling_index=-2.1
   fluence_limit_rshannon = 26
   frb_rate_rshannon = 37
   spectral_index_rshannon = -2.1
   freq_shannon = 1400.00
   point_rshannon = '+'

   # using Ryan Shannon et al, (2018)  spectral index from ASKAP :
   (fluence_rshannon,frb_rates_rshannon) = calc_frb_rates_vs_fluence( options.freq_mhz, fluence_limit_rshannon, frb_rate_rshannon, spectral_index=spectral_index_rshannon, freq_ref=freq_shannon )
   ax_rshannon =  plt.plot( fluence_rshannon,frb_rates_rshannon, linestyle='dashed', color='orange', linewidth=2, markersize=12 )
   plot_list.append(ax_rshannon[0])
   legend_list.append(r'Shannon et al.,2018, $\alpha$=-2.1')
   
   # add Ryan's point iteself (without any frequency scaling, i.e. spectral index = 0 )
   ax_rshannon_point = plt.plot( [fluence_limit_rshannon] , [frb_rate_rshannon] , point_rshannon , color='orange', markersize=20 )
   plot_list.append(ax_rshannon_point[0])
   legend_list.append(r'Shannon et al.,2018, $\alpha$=0')

   # using spectral index = -1
   (fluence_rshannonM1,frb_rates_rshannonM1) = calc_frb_rates_vs_fluence( options.freq_mhz, fluence_limit_rshannon, frb_rate_rshannon, spectral_index=-1 , freq_ref=freq_shannon )
   ax_rshannonM1 =  plt.plot( fluence_rshannonM1,frb_rates_rshannonM1, linestyle='-.', color='orange', linewidth=2, markersize=12 )
   plot_list.append(ax_rshannonM1[0])
   legend_list.append(r'Shannon et al.,2018, $\alpha$=-1')
#   legend_list.append("Shannon et al.,2018, #alpha=-1")

   # using flat (0) spectral index :
   (fluence_rshannon0,frb_rates_rshannon0) = calc_frb_rates_vs_fluence( options.freq_mhz, fluence_limit_rshannon, frb_rate_rshannon, spectral_index=0 , freq_ref=freq_shannon )
   ax_rshannon0 =  plt.plot( fluence_rshannon0,frb_rates_rshannon0, linestyle='dotted', color='orange', linewidth=2, markersize=12 )
   plot_list.append(ax_rshannon0[0])
   legend_list.append(r'Shannon et al.,2018, $\alpha$=0')
#   legend_list.append("Shannon et al.,2018, #alpha=0")
   ##################################################################################################################################################################

   ##################################################################################################################################################################
   #                                     Lovell telescop rate is similar to GBT < 5500 FRBs / day / sky at 332 MHz (not plotted) 
   # /home/msok/Desktop/GRANTS/2021/LIEF/SECTIONS/references/FRB/FRB_Newsletter_202203/2203.04890.pdf 
   ##################################################################################################################################################################

   ##################################################################################################################################################################
   #                                     Parent et al (2020) /home/msok/Desktop/GRANTS/2021/LIEF/SECTIONS/references/FRB/FRB_rates/GBT/acroread/Parent_et_al_2022.pdf
   # Rate = 3400+15400-3300 above fluence 1.44 Jy ms 
   #
   # How to convert flux limit to fluence limit ? Is it assuming Gaussian of peak = Peak flux ?
   #   
   ##################################################################################################################################################################
   fluence_limit_parent_gbt = 1.44
   frb_rate_parent_gbt = 3400
   spectral_index_parent_gbt  = 0. # 350 MHz is very close to we can assume pessimistic (specrtal index = 0) scenario
   freq_parent_gbt = 350.00
   point_parent_gbt = 'x'

   (fluence_parent_gbt,frb_rates_parent_gbt) = calc_frb_rates_vs_fluence( options.freq_mhz, fluence_limit_parent_gbt, frb_rate_parent_gbt, spectral_index=-2, freq_ref=freq_parent_gbt )
   ax_parent_gbtM2 =  plt.plot( fluence_parent_gbt, frb_rates_parent_gbt, linestyle='dotted', color='green', linewidth=2, markersize=12 )
   plot_list.append(ax_parent_gbtM2[0])
   legend_list.append(r'Parent et al.,2020, $\alpha$=-2')

   (fluence_parent_gbt,frb_rates_parent_gbt) = calc_frb_rates_vs_fluence( options.freq_mhz, fluence_limit_parent_gbt, frb_rate_parent_gbt, spectral_index=-1, freq_ref=freq_parent_gbt )
   ax_parent_gbtM1 =  plt.plot( fluence_parent_gbt,frb_rates_parent_gbt, linestyle='dotted', color='green', linewidth=2, markersize=12 )
   plot_list.append(ax_parent_gbtM1[0])
   legend_list.append(r'Parent et al.,2020, $\alpha$=-1')

   (fluence_parent_gbt,frb_rates_parent_gbt) = calc_frb_rates_vs_fluence( options.freq_mhz, fluence_limit_parent_gbt, frb_rate_parent_gbt, spectral_index=0.00, freq_ref=freq_parent_gbt )
   lower_error = numpy.ones(len(fluence_parent_gbt))*3300
   upper_error = numpy.ones(len(fluence_parent_gbt))*15400
   asymmetric_error = [lower_error, upper_error]
   ax_parent_gbt =  plt.plot( fluence_parent_gbt,frb_rates_parent_gbt, linestyle='dotted', color='green', linewidth=2, markersize=12 )
#   plt.errorbar( fluence_parent_gbt, frb_rates_parent_gbt, yerr=asymmetric_error, fmt='.', color='green')
   plt.fill_between( fluence_parent_gbt, frb_rates_parent_gbt-lower_error, frb_rates_parent_gbt+upper_error, color='green' , alpha=0.1 )
   plot_list.append(ax_parent_gbt[0])
   legend_list.append(r'Parent et al.,2020, $\alpha$=0')
   
   # add GBT point itself -  (without any frequency scaling, i.e. spectral index = 0 )
   plt.plot( [fluence_limit_parent_gbt] , [frb_rate_parent_gbt] , point_parent_gbt , color='green', markersize=20 )

   ##################################################################################################################################################################
   
   ##################################################################################################################################################################
   # CHIME :   /home/msok/Desktop/GRANTS/2021/LIEF/SECTIONS/references/FRB/FRB_rates/CHIME 
   #             https://ui.adsabs.harvard.edu/abs/2021ApJ...923....1P/abstract
   # 
   #          In this assumption, the derived FRB rate was RFRB = 820 /sky/day above a fluence F > 5 Jy ms. in /home/msok/Desktop/GRANTS/2021/LIEF/SECTIONS/references/FRB/FRB_Newsletter_202203/2203.04890.pdf 
   #              https://ui.adsabs.harvard.edu/abs/2021ApJS..257...59C/abstract
   # " We infer a sky rate of [820 +/- 60(stat.) +220 -200(sys.)]/sky/day above a fluence of 5 Jy ms at 600 MHz, with a scattering time at 600 MHz under 10 ms and DM above 100 pc cm-3. "
   #
   ##################################################################################################################################################################
   ##################################################################################################################################################################
   fluence_limit_chime = 5
   frb_rate_chime = 820
   freq_chime = 600.00
   point_chime = 'D'
   chime_color='red'

   # 
   (fluence_chime,frb_rates_chime) = calc_frb_rates_vs_fluence( options.freq_mhz, fluence_limit_chime, frb_rate_chime, spectral_index=0, freq_ref=freq_chime )
   ax_chime =  plt.plot( fluence_chime, frb_rates_chime, linestyle='-', color=chime_color, linewidth=2, markersize=12 )
   plot_list.append(ax_chime[0])
   legend_list.append(r'CHIME, 600 MHz, Amiri et al., 2021, $\alpha$=0')
   
   # add Ryan's point iteself (without any frequency scaling, i.e. spectral index = 0 )
   ax_chime_point = plt.plot( [fluence_limit_chime] , [frb_rate_chime] , point_chime , color=chime_color, markersize=10 )
   plot_list.append(ax_chime_point[0])
   legend_list.append(r'CHIME, 600 MHz, Amiri et al., 2021')

   # using spectral index = -1
   (fluence_chimeM1,frb_rates_chimeM1) = calc_frb_rates_vs_fluence( options.freq_mhz, fluence_limit_chime, frb_rate_chime, spectral_index=-1 , freq_ref=freq_chime )
   ax_chimeM1 =  plt.plot( fluence_chimeM1,frb_rates_chimeM1, linestyle='--', color=chime_color, linewidth=2, markersize=12 )
   plot_list.append(ax_chimeM1[0])
   legend_list.append(r'CHIME, 600 MHz, Amiri et al., 2021, $\alpha$=-1')
#   legend_list.append("Shannon et al.,2018, #alpha=-1")

   # using flat (0) spectral index :
   (fluence_chimeM2,frb_rates_chimeM2) = calc_frb_rates_vs_fluence( options.freq_mhz, fluence_limit_chime, frb_rate_chime, spectral_index=-2, freq_ref=freq_chime )
   ax_chimeM2 =  plt.plot( fluence_chimeM2, frb_rates_chimeM2, linestyle='-.', color=chime_color, linewidth=2, markersize=12 )
   plot_list.append(ax_chimeM2[0])
   legend_list.append(r'CHIME, 600 MHz, Amiri et al., 2021, $\alpha$=-2')
   
   lower_error = numpy.ones(len(fluence_chime))*200
   upper_error = numpy.ones(len(fluence_chime))*220
   asymmetric_error = [lower_error, upper_error]
   plt.fill_between( fluence_chime, frb_rates_chime-lower_error, frb_rates_chime+upper_error, color=chime_color , alpha=0.1 )


   ##################################################################################################################################################################
   # LOFAR :   /home/msok/Desktop/GRANTS/2021/LIEF/SECTIONS/references/FRB/FRB_rates/LOFAR/2012.08348_MS.pdf
   #           pm = Pastor-Marazuela, Ines : https://ui.adsabs.harvard.edu/abs/2021Natur.596..505P/abstract
   #           FRB rate based on FRB20180916B burst rate at 150 megahertz, we find there are 3-450 FRBs in the sky per day above 50 Jy ms.
   #           By extending this limit with Euclidean fluence scaling, the rate becomes RFRB = 90  1400 sky1 day1.
   ##################################################################################################################################################################
   fluence_limit_lofar_pm_low = 50
   frb_rate_lofar_pm_low = 3   
   fluence_limit_lofar_pm_high = 50
   frb_rate_lofar_pm_high = 450
   frb_rate_lofar_pm = (frb_rate_lofar_pm_low + frb_rate_lofar_pm_high) / 2.00
   fluence_limit_lofar_pm = (fluence_limit_lofar_pm_low+fluence_limit_lofar_pm_high)/2.00
   
   spectral_index_lofar_pm  = 0. # 150 MHz is very close to we can assume pessimistic (specrtal index = 0) scenario
   freq_lofar_pm = 150.00
   point_lofar_pm = 's'

   (fluence_lofar_pm,frb_rates_lofar_pm) = calc_frb_rates_vs_fluence( options.freq_mhz, fluence_limit_lofar_pm, frb_rate_lofar_pm, spectral_index=0, freq_ref=freq_lofar_pm )
   ax_parent_lofar_pm_sp0 =  plt.plot( fluence_lofar_pm,frb_rates_lofar_pm, linestyle='-', color='blue', linewidth=2, markersize=12 )
   plot_list.append(ax_parent_lofar_pm_sp0[0])
   legend_list.append(r'Pastor-Marazuela, (Lofar) et al., $\alpha$=0')
   
   lower_error = numpy.ones(len(fluence_parent_gbt))*frb_rate_lofar_pm_low
   upper_error = numpy.ones(len(fluence_parent_gbt))*frb_rate_lofar_pm_high
   asymmetric_error = [lower_error, upper_error]
   plt.fill_between( fluence_lofar_pm, lower_error, upper_error,  color='blue', alpha=0.1 )

   ##################################################################################################################################################################
   # LOFAR OLD : /home/msok/Desktop/GRANTS/2021/LIEF/SECTIONS/references/FRB/FRB_Newsletter_202203/2203.04890.pdf
   # 
   #   LPPS : < 150/day/sky at F > 107 Jy ms at frequency ??? MHz ???
   #   Rawlings Array : 29/day/sky at F > 310 Jy ms at 145 MHz 
   #   FRATS : <1400/day/sky F > 6000 Jy ms  between 119 and 151 MHz    
   ##################################################################################################################################################################

   ##################################################################################################################################################################
   # UTMOST : /home/msok/Desktop/GRANTS/2021/LIEF/SECTIONS/references/FRB/FRB_Newsletter_202203/2203.04890.pdf
   #  1/ RFRB < 1000/sky/day at 843 MHz above a fluence F > 23 Jy ms
   #  2/ RFRB = 78 /sky/day above a fluence F > 11 Jy ms at 800 MHz 
   #  3/ RFRB = 98 /sky/day above a fluence F > 8 Jy ms at 800 MHz 
   ##################################################################################################################################################################
   fluence_limit_utmost = 8
   frb_rate_utmost = 98
   freq_utmost = 800.00
   point_utmost = 'P'
   utmost_color='pink'

   # 
   (fluence_utmost,frb_rates_utmost) = calc_frb_rates_vs_fluence( options.freq_mhz, fluence_limit_utmost, frb_rate_utmost, spectral_index=0, freq_ref=freq_utmost )
   ax_utmost =  plt.plot( fluence_utmost, frb_rates_utmost, linestyle='-', color=utmost_color, linewidth=2, markersize=12 )
   plot_list.append(ax_utmost[0])
   legend_list.append(r'UTMOST,???, $\alpha$=0')
   
   # add Ryan's point iteself (without any frequency scaling, i.e. spectral index = 0 )
   ax_utmost_point = plt.plot( [fluence_limit_utmost] , [frb_rate_utmost] , point_utmost , color=utmost_color, markersize=20 )
   plot_list.append(ax_utmost_point[0])
   legend_list.append(r'UTMOST, 800 MHz')

   # using spectral index = -1
   (fluence_utmostM1,frb_rates_utmostM1) = calc_frb_rates_vs_fluence( options.freq_mhz, fluence_limit_utmost, frb_rate_utmost, spectral_index=-1 , freq_ref=freq_utmost )
   ax_utmostM1 =  plt.plot( fluence_utmostM1,frb_rates_utmostM1, linestyle='--', color=utmost_color, linewidth=2, markersize=12 )
   plot_list.append(ax_utmostM1[0])
   legend_list.append(r'UTMOST,???, $\alpha$=-1')
#   legend_list.append("Shannon et al.,2018, #alpha=-1")

   # using flat (0) spectral index :
   (fluence_utmostM2,frb_rates_utmostM2) = calc_frb_rates_vs_fluence( options.freq_mhz, fluence_limit_utmost, frb_rate_utmost, spectral_index=-2, freq_ref=freq_utmost )
   ax_utmostM2 =  plt.plot( fluence_utmostM2, frb_rates_utmostM2, linestyle='-.', color=utmost_color, linewidth=2, markersize=12 )
   plot_list.append(ax_utmostM2[0])
   legend_list.append(r'UTMOST,???, $\alpha$=-2')


#   plt.legend( legend_list, loc=legend_location, fontsize=20)
   plt.legend( plot_list, legend_list, loc=legend_location, fontsize=10)

   # mark region : https://matplotlib.org/3.1.1/api/_as_gen/matplotlib.pyplot.html
   # plt.axhspan(ymin, ymax[, xmin, xmax])   
#   plt.axhspan( 0.0001, 200, 90, 1000, alpha=0.3, facecolor='0.5' )
   plt.axvspan( 200, 10000, ymin=0.235, ymax=0.535, alpha=0.5, facecolor='red' )
   # legend of the span :
   plt.axvspan( 10, 600, ymin=0.94, ymax=0.99, alpha=0.5, facecolor='red' )
   plt.text(15,14000000,"FRB rate below 350 MHz based on the literature",fontsize=15)


   
   plt.xlabel('Fluence [Jy ms]' , fontsize=30 )
   plt.ylabel('#FRBs / day / sky' , fontsize=30 )
   plt.xscale('log')
   plt.yscale('log')

   plt.savefig( "frb_rates.png" )
   plt.show()   



if __name__ == "__main__":
    (options, args) = parse_options(1)
    
    do_plots(options)
    
        
