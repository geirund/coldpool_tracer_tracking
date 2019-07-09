import numpy as np
import numpy.ma as ma
import subprocess
import os   # for environment 
import math
import matplotlib
import operator
import matplotlib.pyplot as plt
import datetime
from operator import itemgetter
from netCDF4 import Dataset
from six.moves import cPickle   # to save class files
import argparse

atime = datetime.datetime.now()
###############################################################
# SETTINGS
# ###############################################################
DIR = '/scratch/snx3000/geirund/celltrack/celltrack_olga/'
#EXPID='test1plus4K_finneu1/'
parser = argparse.ArgumentParser(description=('open file and do whatever'))
parser.add_argument('filename', type=str,help=('file name'))
args = parser.parse_args()
EXPID=args.filename

# READ HEADER FILE TO GET SPLITTING
f = open(DIR+EXPID+'/output/raincell/headerfile.txt', 'r')
lines = f.readlines()
f.close()
merger={}
dur={}
trackStart={}
for line in lines:
  columns = line.split()
  ID = int(columns[0])              #Track ID
  durat=(int(columns[2]))           #duration of track
  starting = (int(columns[3]))      #situation at track beginning 
#  if starting != 0:                # 0 emerging by itself
#   merger[ID] =0
  dur[ID] =  durat
  trackStart[ID] = int(columns[1])

# READ TRACK FILE AND GET TERMINATING OBJ ID
f = open(DIR+EXPID+'/output/raincell/irt_objects_output.txt', 'r')
lines = f.readlines()
f.close()
next_objID={}
for line in lines:
  columns = line.split()
  objID_in      = int(float(columns[2]))
  next_objID_in = int(float(columns[18]))  # object of this track at the next timestep
  next_objID[objID_in] = next_objID_in
#  print 'for obj',objID_in, 'nextobj', next_objID_in 


# READ TRACK FILE AND GET TERMINATING OBJ ID
f = open(DIR+EXPID+'/output/raincell/irt_tracks_output_pure.txt', 'r')
lines = f.readlines()
f.close()
TrackID_forObj = {}
startTracer={}
endObjp1={}
ta = 0
for line in lines:
  columns = line.split()
  ID = int(float(columns[0]))    # Track ID
  ta = max(ID,ta)
  objID_in   = int(float(columns[3]))  # OBJ ID 
  time_in    = int(columns[1])      # time
  conv_in     = float(columns[6])    # divergence
  size_in     = float(columns[4])    # size
  print(conv_in)
  age = time_in-trackStart[ID] +1 
#  if conv_in < -0.0025 and size_in > 50 and not ID in startTracer.keys():
     # first time of this track, wenn divergence > 0.025 ->start tracer
  startTracer[ID] = time_in
#     print conv_in
  TrackID_forObj[objID_in]=ID  # trrck ID for all objects beginning a new track

  # last object of this track
  if age == dur[ID]:
    print(ID, objID_in)
    endObjp1[ID] = next_objID[objID_in] # OBJ ID of next timestep, which is zero if terminates or the ID within the trackit merges into
#    if not ID in startTracer.keys(): # if divergence never high enough, dont set tracer, set merging to 0
#      merger[ID] = 0
#      startTracer[ID] = 500


for ID in endObjp1.keys():
    if ID not in merger.keys():
      if endObjp1[ID] == 0: # no following object to this track -> not merging into another
        merger[ID] = ID
      else:
        if endObjp1[ID] not in TrackID_forObj.keys():
          merger[ID] = 0
        else:
          merger[ID] =  TrackID_forObj[endObjp1[ID] ]  # t




f = open(DIR+EXPID+'/input/cp/mergingCPs.txt','w+')
 
for ID in sorted(merger.keys()):
#  print ID, merger[ID], trackStart[ID], startTracer[ID] 
# just start tacer at start of precip at the moment
  value =str(str(ID) + ' ' +str(merger[ID]) +' ' + str(trackStart[ID])+' '+ str(trackStart[ID]))

#  value =str(str(ID) + ' ' +str(merger[ID]) +' ' + str(startTracer[ID])+' '+ str(trackStart[ID]))
  f.write("%s\n" % value)

os.system("echo "+str(ta) +" > "+DIR+EXPID+"/input/cp/na.txt")

