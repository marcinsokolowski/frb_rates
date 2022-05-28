
from optparse import OptionParser
import math
import sys

def parse_options(idx):
   global debug_level

   usage="Usage: %prog [options]\n"
   usage += "Memory requirements for a specified dispersion delay :\n"
   
   parser = OptionParser(usage=usage,version=1.00)

   # TEST : azim=0, za=0, lst=15.4 
   parser.add_option('-i','--inttime','--integration_time',dest="integration_time",default=0.01, help="Integration time in seconds [default %default seconds]",type="float")
   parser.add_option('-b','--bw','--bandwidth',dest="bandwidth_mhz",default=31.25, help="Observing bandwidth in MHz [default %default MHz]",type="float")
   parser.add_option('-s','--sefd','--telescope_sefd','--station_sefd',dest="sefd",default=2500, help="System Equivalent Flux Density (SEFD) in Jy [default %default Jy]",type="float")
   parser.add_option('-p','--n_pols','--npols','--np','--number_polarisations',dest="npol",default=2, help="Number of polarisations [default %default]",type="int")
   parser.add_option('-c','--n_chan','--nchan','--nc','--number_channels','--n_channels',dest="nchan",default=40, help="Number of channels [default %default]",type="int")
   parser.add_option('--size','--pixels','--px','--npx',dest="n_pixels",default=200, help="Image size in pixels [default %default]",type="int")
   parser.add_option('-d','--dispdelay','--dispersion_delay',dest="dispersion_delay",default=38.9, help="Dipersion delay in seconds [default %default seconds]",type="float")
   parser.add_option('--pixel_size','--n_bytes','--nbytes','--bytes_per_pixel',dest="bytes_per_pixel",default=4, help="Number of bytes per pixel [default %default]",type="int")
   # parser.add_option('-T','--tingay','--show_tingay',action="store_true",dest="show_tingay",default=False, help="Show Tingay et al. (2015) upper limit [default %default]")
   # parser.add_option('-f','--freq','--freq_mhz',dest="freq_mhz",default=200.00, help="Survey frequency in MHz [default %default MHz]",metavar="float",type="float")
   # parser.add_option('--include_alpha_minus2','--alpha_minus2',action="store_true",dest="include_alpha_minus2",default=False, help="Include curves for spectral index alpha -2 [default %default]")
   (options, args) = parser.parse_args(sys.argv[idx:])

   print("#####################################################################################################")
   print("PARAMETERS:")
   print("#####################################################################################################")
   print("Integration time = %.6f [sec] = %.3f [ms]" % (options.integration_time,options.integration_time*1000.00))
   print("Bandwidth = %.6f [MHz] = %.2f [Hz]" % (options.bandwidth_mhz,options.bandwidth_mhz*1000000.00))
   print("SEFD = %.6f [Jy]" % (options.sefd))
   print("Number of polarisations = %d" % (options.npol))
   print("Dipersion delay = %.4f [sec]" % (options.dispersion_delay))
   print("Image size = %d x %d" % (options.n_pixels,options.n_pixels))
   print("Number of bytes per pixel = %d" % (options.bytes_per_pixel))
   print("#####################################################################################################")

   return (options, args)


def check_memory( options ):
   npol = options.npol
   integration_time = options.integration_time
   integration_time_ms = integration_time*1000.00
   bandwidth_hz = options.bandwidth_mhz*1000000.00
   n_channels = options.nchan

   n_images_per_sec = (1.00/options.integration_time)*n_channels
   image_size = options.n_pixels*options.n_pixels*options.bytes_per_pixel # in bytes 
   
   memsize_bytes = image_size*n_images_per_sec*options.dispersion_delay
   memsize_mb = float(memsize_bytes)/1000000.00
   memsize_gb = float(memsize_bytes)/1000000000.00
   mem_bw = memsize_gb/options.integration_time
   
   print("Required memory size = %d bytes = %.3f MB = %.6f GB" % (memsize_bytes,memsize_mb,memsize_gb))
   print("Required memory bandwidth = %.3f GB/s" % (mem_bw))
      
   
   
if __name__ == "__main__":
   (options, args) = parse_options(1)
   
   check_memory( options )   
   


