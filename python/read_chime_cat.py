import pandas as pd
import xlrd 
import csv
# import json
import cPickle as pickle

import numpy
import sys
import matplotlib.pyplot as pyplot
import copy
import math
import os,re,errno
from lmfit import *
from lmfit.models import GaussianModel,LinearModel,PolynomialModel
from optparse import OptionParser,OptionGroup


def parse_options():
   usage="Usage: %prog [options]\n"
   usage+='\tRead and plot CHIME catalogue data\n'
   parser = OptionParser(usage=usage,version=1.00)
   parser.add_option('-d','--device_id','--device',dest="device_id",default=None, help="Device ID [default %default]")
   parser.add_option('-t','-c','--crop_type',dest="crop_type",default="Wheat", help="Plant type [default %default]")
   parser.add_option('-o','--outfile','--outf',dest="output_file",default="results.txt", help="Results outputfile [default %default]")
   parser.add_option('--filecolumn',dest="filecolumn",default=" ", help="Column with file name in csv file [default %default]")
   parser.add_option('-i','--infile','--inf',dest="input_file",default="N16WC1_TOS1.csv", help="Input CSV file [default %default]")
   parser.add_option('--outdir',dest="outdir",default=None, help="Output dir [default %default]")
   parser.add_option('-f','--dofit','--fit',dest="dofit",default=False,action="store_true", help="Perform fitting [default %default]")
   parser.add_option('-p','--doplot','--plot',dest="doplot",default=False,action="store_true", help="Do plot [default %default]")
   parser.add_option('-x','--scan_number','--plant_index','--number',dest="scan_number",default=None,help="Scan number, example 1,3 or just 1 [default %default]") # example "1,3"
   parser.add_option('-b','--barcode_index','--barindex','--index',dest="barcode_index",default=None,help="Barcode index [default %default]",type="int")
   parser.add_option('--offset2channel','--offset','--offset_channel2value',dest="offset_channel2value",default=None,help="Offset data so that channel 0 has this value [default %default]",type="float")
   parser.add_option('--cut_n_sigma','--n_sigma_cut',dest="n_sigma_cut",default=None,help="Cut outlier spectra more then N sigma away from mean [default %default]",type="float")
   parser.add_option('-r','--refl','--reflectance',dest="reflectance",default=False,action="store_true", help="Reflectance [default %default]")
   parser.add_option('--derv','--derivative',dest="derivative",default=False,action="store_true", help="Calculate 1st derivative [default %default]")
   parser.add_option('-s','--site_name',dest="site_name",default=None, help="Site name [default %default]")   
   parser.add_option('-a','--average','--avg',dest="do_average",default=None,help="Average data [default %default]")
   
   parser.add_option('--prefix',dest="prefix",default="DL",help="Requested prefix [default %default]")
   parser.add_option('--pga_gain',dest="pga_gain",default=32,help="PGA gain [default %default]",type="int")

   parser.add_option('--mean_spectrum',dest="mean_spectrum",default=None, help="Mean spectrum [default %default]")
   parser.add_option('--rms_spectrum',dest="rms_spectrum",default=None, help="RMS spectrum [default %default]")

   parser.add_option('--data_format',dest="data_format",default="dry_fresh_lab", help="Data format [default %default]")
   parser.add_option('--excel2csv','--convert',dest="excel2csv",default=None, help="Convert excel file to csv [default %default]")
   

#   parser.add_option('-p','--pol',dest="selected_pol",default="X", help="Selected polarisation",metavar="STRING")
#   parser.add_option('-f','--freq',dest="freq",default=230, help="Frequency in MHz",metavar="FLOAT",type="float")
#   parser.add_option('-d','--debug',action="store_true",dest="debug",default=False, help="Debug mode")
#   parser.add_option('-v','--view',action="store_true",dest="view_only",default=False, help="View only (no saving)")
   (options, args) = parser.parse_args()
   return (options, args)   

def print_options(options,args) :

   scan_list = None 
   if options.scan_number is not None :
      scan_list = options.scan_number.split(",")
      scan_list = map(int,scan_list)

   print "####################################################"
   print "PARAMTERS :"
   print "####################################################"
   print "####################################################"

   
   return (options, args)



def mkdir_p(path):
   try:
      os.makedirs(path)
   except OSError as exc: # Python >2.5
      if exc.errno == errno.EEXIST:
         pass
      else: raise
       
       
def is_float(input):
  try:
    num = float(input)
  except ValueError:
    return False
  return True

def is_number(s):
    """ Returns True is string is a number. """
    try:
        int(s)
        return True
    except ValueError:
        return False



def read_spectrum( file_name ) :

   data_x=[]
   data_y=[]
   file=open( file_name,'r')
   
   # reads the entire file into a list of strings variable data :
   data=file.readlines()
   for line in data : 
       words = line.split(' ')

       if line[0] == '#' :
           continue

       if line[0] != "#" :
           x=float(words[0+0])
           y=float(words[1+0])
      
       data_x.append(x)
       data_y.append(y)

       file.close()
   
   print "Read %d points from file %s" % (len(data_x),file_name)
   return (data_x,data_y)    


def excel2csv( xls_file = "test.xlsx" , csv_file="test.csv", sheet_name="test") :
#    file="140818_Combined_Data_for_Curtin_2018.csv"
    wb = xlrd.open_workbook( xls_file )
    sheet = wb.sheet_by_name( sheet_name )

    csv_f = open(csv_file,"wb")
    wr = csv.writer( csv_f, quoting=csv.QUOTE_ALL)
    for rownum in xrange(sheet.nrows) :
        wr.writerow(sheet.row_values(rownum))
    csv_f.close()




if __name__ == '__main__':
   (options, args) = parse_options()
   print_options( options, args) 

   file = "chimefrbcat1.csv"
   df = pd.read_csv(file,low_memory=False)

   for i in range(0,df.shape[0]) :
      print("%s %.4f %.4f %.4f" % (df['tns_name'][i],df['dm_fitb'][i],df['fluence'][i],df['fluence_err'][i]))
      