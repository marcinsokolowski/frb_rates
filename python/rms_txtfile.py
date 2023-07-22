# based on examples
# http://docs.astropy.org/en/stable/io/fits/ 
# https://python4astronomers.github.io/astropy/fits.html

import astropy.io.fits as pyfits
import pylab
import math
import numpy
from array import *
import matplotlib.pyplot as plt
import numpy as np
import string
import sys
import os
import errno
import getopt
from optparse import OptionParser,OptionGroup
# from histofile import histoarray


def parse_options():
   usage="Usage: %prog [options]\n"
   usage+='\tCalculate MEAN and RMS of values in text file\n'
   parser = OptionParser(usage=usage,version=1.00)
#   parser.add_option('-t','--threshold',dest="threshold",default=5, help="Threshold to find sources (expressed in sigmas) [default %default sigma]",type="float")
#   parser.add_option('-f','--fits_type',dest="fits_type",default=0, help="FITS file type 0-single x,y image, 1-multi-axis [default %default]",type="int")
#   parser.add_option('-w','--window',nargs = 4, dest="window", help="Window to calculate RMS [default %default]",type="int")   
#   parser.add_option('-r','--regfile',action="store_true",dest="save_regfile",default=False, help="Save regfile with found pixels [default %default]")
#   parser.add_option('-x','--x',dest="position_x",default=-1,help="X coordinate to calculate RMS around [default %default]",type="int")
#   parser.add_option('-y','--y',dest="position_y",default=-1,help="Y coordinate to calculate RMS around [default %default]",type="int")
#   parser.add_option('--radius',dest="radius",default=1,help="Radius to calculate RMS around (X,Y) passed in -x and -y options [default %default]",type="int")
#   parser.add_option('--plotname',dest="plotname",default=None,help="Png file name if plot is required too [default %default]",type="string")
   (options, args) = parser.parse_args()
   return (options, args)



if __name__ == '__main__':
   filename='data.txt'
   if len(sys.argv) > 1:
      filename = sys.argv[1]
       
   file=open(filename,'r')
   
   debug = False
   

   # reads the entire file into a list of strings variable data :
   data=file.readlines()
   # print data

   # initialisation of empty lists :
   x_arr=[]
   y_arr=[]
   # x_err=[]
   # y_err=[]

   test_x=0
   test_y=0
   count=0

   for line in data : 
      words = line.split(' ')

      if line[0] == '#' :
         continue

      if line[0] != "#" :
         x=float(words[0+0])
         y=float(words[1+0])
      
         x_arr.append(x)
         y_arr.append(y)
      
         test_x += x
         test_y += y 
         count += 1
      
   mean_x = numpy.mean(x_arr)
   mean_y = numpy.mean(y_arr)
   rms_x  = numpy.std(x_arr)
   rms_y  = numpy.std(y_arr)

   if debug : 
      print "------------------------------------------------------------------ MEAN  ------------------------------------------------------------------"
      print "mean_x = %020.10f (TEST = %020.10f)" % (mean_x,test_x/count)
      print "mean_y = %020.10f (TEST = %020.10f)" % (mean_y,test_y/count)
   
      print "------------------------------------------------------------------ RMS  ------------------------------------------------------------------"
      print "rms_x = %020.10f" % (rms_x)
      print "rms_y = %020.10f" % (rms_y)



   # print "------------------------------------------------------------------ VAR  ------------------------------------------------------------------"
   # print "var_x = %020.10f" % (numpy.var(x_arr))
   # print "var_y = %020.10f" % (numpy.var(y_arr))

   # print "------------------------------------------------------------------ PEARSONR  ------------------------------------------------------------------"
   # print "pearsonr(x,y) = %020.10f" % (scipy.stats.pearsonr(x,y))
   # version 1:
   # corr =  pearsonr(x_arr,y_arr)
   # version 2 :
   # corr = stats.pearsonr(x_arr,y_arr)
   # print corr

   # print "------------------------------------------------------------------ LINEGRESS  ------------------------------------------------------------------"
   # version 1 :
   # slope,intercept, rval, pval, stderr = linregress(x_arr,y_arr)
   # version 2 :
   # slope,intercept, rval, pval, stderr = stats.linregress(x_arr,y_arr)
   # print "y = {0}*x + {1}".format(slope,intercept)
   # print "rval = {0}".format(rval)
   # print "pval = {0}".format(pval)
   # print "stderr = {0}".format(stderr)
         
   print("FRBs per day  : %.4f +/- %.4f" % (mean_y,rms_y))
   print("FRBs per year : %.4f +/- %.4f" % (mean_y*365,rms_y*365))
   