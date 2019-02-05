import numpy as np
import numpy.ma as ma
import os   # for environment 
import math
import matplotlib
import operator
import matplotlib.pyplot as plt
import datetime
from classes import TRACER,CP_map,  CPstart, CPlife
from operator import itemgetter 
from netCDF4 import Dataset   
from six.moves import cPickle   # to save class files

atime = datetime.datetime.now()

###############################################################
# SETTINGS
# ###############################################################
#EXPID ='totest' #multismall'
EXPID ='lindp2K_fin4' #'lindp2K_fin' #test1plus4K_fin3' #lindp2K_fin' #

lcol = False   # no collison information if True
lcpstart = True
lmap = False
lmap2 = True #faster

odir = '/nbi/ac/conv1/henneb/results/coldpool/'

############################################
# INITIALIZE FIELDS AND READ DATA
#############################################
data={}
xylist={}   # list to check if map at x|y in key time needs to be initialized

f = open(odir+EXPID+'/output/cp/numberoftracer.txt', 'r')
lines = f.readlines()
column=lines[0].split()
nooftracer = int(column[0])+1
for i in range(1,nooftracer):
  data[i] = TRACER()

# timing analysis:
timedummy1 = datetime.datetime.now()
timedummy2 = datetime.datetime.now()
sumtimet = {}
nop = {}

if lcpstart:
    CPinit = {}
    timedummy1 = datetime.datetime.now()
    timedummy2 = datetime.datetime.now()
    sumlcpstarttime = timedummy2 - timedummy1

    #for line in lines:
    #  columns = line.split()
    #  raincell= (int(columns[0])) 
    #  rainmerger= (int(columns[1]))
    #  if rainmerger != 0:
    #    print raincell
    #    data[raincell] = TRACER()

if lmap or lmap2:
    CoPoMap = {}
    timedummy1 = datetime.datetime.now()
    timedummy2 = datetime.datetime.now()
    sumlmaptime = timedummy2 - timedummy1
    xylist = {}

column=lines[0].split()
timeold = int(column[0])

timeold = int(column[0])
sumtimet[timeold] = timedummy2-timedummy1
nop[timeold] = 0
xylist[timeold] = {} #[]

f = open(odir+EXPID+'/output/cp/coldpool_tracer_out_all.txt', 'r')
lines = f.readlines()

##########################################
# START MAIN LOOP
########################################
for line in lines:
    time1 = datetime.datetime.now()
    columns = line.split()
    tist = (int(columns[0]))     # timestep
    if tist != timeold:
        if lmap:
          xylist[tist] =[]
        elif lmap2:
          xylist[tist] = {}
#  if tist < 90: # just for testing and not running trough all data
    age  = (int(columns[1]))     # age
    tID  = (int(columns[2]))     # tracer ID
    cCP  = (int(columns[3]))     # belonging to CP ID
    xpos1= (float(columns[4]))
    ypos1= (float(columns[5]))
    xpos = (int(columns[6]))     # xpos of tracer
    ypos = (int(columns[7]))     # ypos ...
    dist = (float(columns[8]))     # distance to cog
    phi  = (float(columns[9]))     # angle to y achsis
    uvel = (float(columns[10]))    # velocity in x direction
    vvel = (float(columns[11]))    # velocity in y direction
    xd   = (float(columns[12])) 
    yd   = (float(columns[13]))
    cogx = (float(columns[14]))
    cogy = (float(columns[15]))
    mm   = (int(columns[16]))    # merger = 1
    pID  = (int(columns[17]))    # belonging to precip ID (can differ from CP ID for events that merge into another
    #print 'at time', tist ,'for CP', cCP

    FF   = (uvel**2+vvel**2)**0.5
    if uvel > 0 and vvel > 0:
      alpha = math.atan(vvel/uvel)
    if uvel < 0 and vvel >0:
      alpha = math.pi/2. * math.atan(abs(uvel/vvel)) +math.pi/2.
    if uvel < 0 and vvel < 0:
      alpha = math.atan(vvel/uvel) + math.pi
    if uvel > 0 and vvel < 0:
      alpha =  math.atan(uvel/abs(vvel)) + 3./2. *math.pi
    if uvel == 0 and vvel >0:
       alpha = math.pi/2.
    if uvel == 0 and vvel < 0:
       alpha = 3./2.*math.pi
    if vvel == 0 and uvel > 0:
       alpha = 0.
    if vvel == 0 and uvel<0:
       alpha = math.pi

    vr   = uvel *math.cos(phi) + vvel * math.sin(phi)
    vt   = vvel *math.cos(phi) - uvel * math.sin(phi)
    bstoretime = datetime.datetime.now()
    # STORE THE TRACER DATA sorted by tracer ID
    data[tID].add(tist,age,cCP,xpos,ypos,dist,phi,uvel,vvel,FF,vr,vt,xd,yd,cogx,cogy,alpha,xpos1,ypos1,mm,pID)
    astoretime = datetime.datetime.now()
    if lcpstart:
        blcptime = datetime.datetime.now()
        if pID not in CPinit.keys():
           CPinit[pID]=CPstart(tist,cogx,cogy,dist)
        if tist in range(CPinit[pID].time,CPinit[pID].time+3):
          CPinit[pID].r = max(CPinit[pID].r,dist)
        alcpstarttime = datetime.datetime.now()
        sumlcpstarttime += (alcpstarttime-blcptime)

###############################################################################################
# lmap: write information on map  format
##############################################################################################
    if lmap2: #check if this is faster
        blmaptime = datetime.datetime.now()
        if xpos not in xylist[tist].keys():
            CoPoMap[tist,xpos,ypos]=CP_map()
            xylist[tist][xpos] = []
        elif ypos not in xylist[tist][xpos]:
            CoPoMap[tist,xpos,ypos]=CP_map()
            xylist[tist][xpos].append(ypos)
        CoPoMap[tist,xpos,ypos].add(pID,age,tID)
    alcoltime = datetime.datetime.now()
    time2 = datetime.datetime.now()
    if tist == timeold:
        sumtimet[tist]  += time2-time1
        nop[tist] +=1
    else:
        sumtimet[tist] = timedummy2-timedummy1
        print 'at timestep', tist,sumtimet[tist] 

        nop[tist] = 0

    timeold=tist

b = {}
for (t,x,y) in CoPoMap.keys() :
    #print x,y,t     
    cCP = CoPoMap[t,x,y].CPs
 
    # summerize CP lifecycle
    for cp in list(cCP.values()):
      if cp not in b.keys():
        b[cp] = CPlife(cp,t)
      if len(list(cCP.values())) ==1 :
        b[cp].add(CoPoMap[t,x,y].nTrCP[cp],CoPoMap[t,x,y].nTrCP[cp],1,t,CoPoMap[t,x,y].age[cp],x,y)
      else:
        b[cp].add(CoPoMap[t,x,y].nTrCP[cp],0,0,t,CoPoMap[t,x,y].age[cp],x,y)
        for cps in list(cCP.values()):
         if not cps == cp:
          b[cp].add_others(t,cps,CoPoMap[t,x,y].nTrCP[cp],CoPoMap[t,x,y].nTrCP[cps] )
 

###############################################################
# SAVE DATA
############################################################
btime = datetime.datetime.now()

print 'took ',(btime-atime)

#if not os.path.exists(odir+EXPID+'/tempdata/'):
# os.makedirs(odir+EXPID+'/tempdata/')
f = open(odir+EXPID+'/output/cp/Tracer.save', 'wb')
cPickle.dump(data, f, protocol=cPickle.HIGHEST_PROTOCOL)
f.close()

if lmap or lmap2:
  f = open(odir+EXPID+'/output/cp/TracerMap.save', 'wb')
  cPickle.dump(CoPoMap, f, protocol=cPickle.HIGHEST_PROTOCOL)
  f.close()
  print 'lmap took:', sumlmaptime

if lcpstart:
  f = open(odir+EXPID+'/output/cp/CPstart.save', 'wb')
  cPickle.dump(CPinit, f, protocol=cPickle.HIGHEST_PROTOCOL)
  f.close()
  print 'lcpstart took:', sumlcpstarttime

f = open(odir+EXPID+'/output/cp/CPlife.save', 'wb')
cPickle.dump(b, f, protocol=cPickle.HIGHEST_PROTOCOL)
f.close()


etime = datetime.datetime.now()
#for k in sumtimet.keys():
#    print 'for', k,'it took:', sumtimet[k], 'for ',nop[k],'operations', sumtimet[k]/nop[k]



print 'took ',(etime-atime)


