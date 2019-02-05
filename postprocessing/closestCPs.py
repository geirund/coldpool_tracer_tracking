import numpy as np
import numpy.ma as ma
import os   # for environment
import math
import matplotlib
import operator
import copy
import matplotlib.pyplot as plt
import datetime
from classes import TRACER,CP_map,  CPstart, RAINCLOUD, CPlife
from operator import itemgetter
from netCDF4 import Dataset
from six.moves import cPickle   # to save class files
starttime = datetime.datetime.now()
EXPID ='test1plus4K_fin4' #'lindp2K_fin' #'test1plus4K_fin2'
#EXPID ='lindp2K_neworder'
DIR = '/nbi/ac/conv1/henneb/results/test1plus4K/'# lind_p2K/' #test1plus4K/'
# SETTINGS
dir = '/nbi/ac/conv1/henneb/results/coldpool/'

f = open(dir+EXPID+'/output/cp/CPstart.save', 'rb')
CPinit = cPickle.load(f)
f.close()

f = open(dir+EXPID+'/output/cp/TracerMap.save', 'rb')
CPmap = cPickle.load(f)
f.close()

f = open(DIR+'/r_int1/output/timedelay.save', 'rb')
timedel= cPickle.load(f)
f.close()

f = open(dir+EXPID+'/output/cp/Tracer.save', 'rb')
Tracer = cPickle.load(f)
f.close()

f = open(dir+EXPID+'/output/cp/CPlife.save', 'rb')
CPinfos = cPickle.load(f)
f.close()

opentime = datetime.datetime.now()
print 'time to load files:', opentime-starttime

parentCP = {}  # list of parent (causing CP) for kids(caused) event in key
kidsCP = {}    # list of kids (caused CP) by parent(causing) CP in key
class CPparent():
    def __init__(self):
      self.x = {}
      self.y = {}
      self.dist = {}
      self.noGP = {}
      self.age = {}
    def addCP(self,ID,dist):
      self.dist[ID] =  dist
    def delete(self,ID):
      del self.dist[ID]
      if not self.dist.keys():
          del self
    def addinfo(self,ID,noGP,age):
      self.noGP[ID] = noGP
      self.age[ID] = age
#xylist={}
#for (t,x,y) in CPmap.keys():
#   if not t in xylist.keys():
#       xylist[t] = []
#   xylist[t].append((x,y))

txydic = {}
for (t,x,y) in CPmap.keys():
    if not t in txydic.keys():
        txydic[t] = {}
    if not x in txydic[t].keys():
        txydic[t][x]=[]
    txydic[t][x].append(y)
#xylist = {a:(b,c) for a,b,c in CPmap.keys()}
#txydic = {a:{b:c} for a,b,c in CPmap.keys()}
timedummy1 = datetime.datetime.now()
timedummy2 = datetime.datetime.now()
sumif1 = timedummy2 - timedummy1
sumif2 = timedummy2 - timedummy1
sumif3 = timedummy2 - timedummy1
# start main routine
print 'start main loop across all newly formed CPs sorted by time'
timearray= {}
#del CPinit[0]
for i in CPinit.keys():
 if i in timedel.keys():
    dt = CPinit[i].t - timedel[i].delay
    timearray[dt] = CPinit[i].t - timedel[i].delay
for k in set(sorted(list(CPinit[i].t for i in sorted(CPinit.keys())))):
#print timearray.keys()
#for i in range [0,1]:
#for sorted(list(timearray.keys())):
 print 'timestep',k
 for ID in [kk for kk in CPinit.keys() if CPinit[kk].t ==k]: # and CPinit[kk].t > 154 and CPinit[kk].t < 160]:
     if not ID in parentCP.keys():  #make dictonary entry for every CP
         parentCP[ID] = CPparent()
     # consider timedelay from initiation to coldpool formation 
     t = CPinit[ID].t - timedel[ID].delay
     print ID, CPinit[ID].t
     #define a search radius in which cold pools can have affected the freshly formed CP
     xmax = int(1.5 * CPinit[ID].r)
     ymax = int(1.5 * CPinit[ID].r)
     rint = int(CPinit[ID].r)
     print 'init radius', rint
     # search in initial area of CP ring, (reduce search range to make program faster)
     for  i in [ii for ii in range(-1*xmax-1,xmax+1,1)]: #if ii not in range(int(-0.5*rint+1),int(0.5*rint-1),1)]:
      for j in [jj for jj in range(-1*ymax-1,ymax+1,1) if (i**2 + jj**2)**0.5<= rint+1]: # if jj not in range(int(-0.5*rint+1),int(0.5*rint-1),1) and (i**2 + jj**2)**0.5<= rint+1]:
         # gripoint in map to check for tracer existence:
         x = (((int(CPinit[ID].x)+i)+320-1)%320)+1
         y = (((int(CPinit[ID].y)+j)+320-1)%320)+1


         if t in txydic.keys():
           if x in txydic[t].keys():
             if y in txydic[t][x]:
#         if (t,x,y) in CPmap.keys():
             # loop trough all CPs with active tracers at xy but tracer of freshly formed CP (ID)
              for CPs in [cc for cc in CPmap[t,x,y].CPs.values() if cc != ID]:
               if 1 in list(Tracer[it].active[t] for it in list(CPmap[t,x,y].tracers[CPs])):
                 # create dictonary entry for ID
                 if not ID in parentCP.keys():
                   parentCP[ID] = CPparent()
                 if not CPs in kidsCP.keys():
                   kidsCP[CPs] = CPparent()
                 dist = (((CPinit[ID].x- CPinit[CPs].x)**2. + (CPinit[ID].y- CPinit[CPs].y)**2.)**0.5)%320
                 #if CPs has caused one of the others which now also cause ID it can not act as parent event but is a grani.
                 # not a sufficient condition because not at all gps where CP conributes to ID all kids are present,
                 # thus ID might be a valid triggerer at other location
                 if (set(kidsCP[CPs].dist.keys()) & set(CPmap[t,x,y].CPs.values())):
                     print CPs,'has caused', list(set(kidsCP[CPs].dist.keys()).intersection(set(CPmap[t,x,y].CPs.values())))
                 else:
                  parentCP[ID].addCP(CPs,dist)
                  for cpi in parentCP[ID].dist.keys():
                    parentCP[ID].addinfo(cpi,CPinfos[cpi].noGP[t],CPinfos[cpi].age[t])
                  kidsCP[CPs].addCP(ID,dist)
                  print ID, 'was caused by', CPs
                 # set collided tracer inactive:
                 # all at gridbox, x,y belonging to causing CP for all future timesteps
                 for tr in CPmap[t,x,y].tracers[CPs]:
                     for tk in [tkeys for tkeys in Tracer[tr].active.keys() if tkeys >= t]:
                       Tracer[tr].active[tkeys] = 0

#
#delete function:
# delete all CPs causing an event together with an CP they have caused theirself:
#*** to DO : also set according tracers inactive!!!
parentCP2 = copy.deepcopy(parentCP)
kidsCP2 = copy.deepcopy(kidsCP)
parentCP3 = copy.deepcopy(parentCP)
kidsCP3 = copy.deepcopy(kidsCP)
print 'before change'
for k in sorted(parentCP.keys()):
    if k in parentCP2.keys():

        print k, 'is caused by:', parentCP[k].dist.keys()

#print kidsCP[22].dist.keys()
for k in parentCP.keys():  #alle CPs
  for parent in parentCP[k].dist.keys(): 
    if parentCP[k].noGP[parent] < 20:
      parentCP3[k].delete(parent)
# if k == 105:
  if len(list(parentCP[k].dist.keys())) > 1:  # only if more than one CP is included
    for parent in parentCP[k].dist.keys():
        print 'loop', parent, 'in', parentCP[k].dist.keys()
        print parent, 'is mami of:', kidsCP[parent].dist.keys()
        print 'is one of  the kids from ',parent, 'also Mami of current baby?'
        for kid in kidsCP[parent].dist.keys():
            print 'has ',kid, 'also caused the current baby?'
            if kid in parentCP[k].dist.keys():
             print 'yes, Mami', parent, 'is Omi and should be deleted'
            #if kid is causing same event as its parent, parent ate granis and cant cause this event anymore
             #if k in parentCP2.keys():
             print k
             if parent in parentCP2[k].dist.keys():
              # if hasattr(parentCP2,'dist'): # in parentCP2.dist.keys():
               print 'delete parent', parent
               parentCP2[k].delete(parent)
                 #delete parent[] if all kids were removed . can this happen?
                 #if not parentCP2[k].delete
             #if parent in kidsCP2.keys():
             #  if hasattr(kidsCP2,'dist'): # in kidsCP2.dist.keys():
             if k in kidsCP2[parent].dist.keys():

              print 'delete kid'
              kidsCP2[parent].delete(k), k

###############################################
# SAVE DATA
################################################
# read merger list
f = open('/nbi/ac/conv1/henneb/results/coldpool/test1plus4K/input/cp/mergingCPs.txt', 'r')
lines = f.readlines()
merger = {}
for line in lines:
    columns = line.split()
    rID = (int(columns[0])) # original ID
    mID = (int(columns[1])) # merger ID
    merger[rID] = mID
for k in sorted(parentCP.keys()):
    print 'no of gp', parentCP[k].noGP.values()
    print 'age of gp', parentCP[k].age.values()    
    print k, 'is caused by 1:', parentCP[k].dist.keys()
    print k, 'is caused by 2:', parentCP2[k].dist.keys()  , 'generation rule'
    print k, 'is caused by 3:', parentCP3[k].dist.keys()  ,' less than 10 gps'

#     if k in parentCP2.keys():
#        print k, 'is caused by:', parentCP2[k].dist.keys(), parentCP[k].dist.keys()
#    else:
#        print k, 'is caused by:', '[ ]', parentCP[k].dist.keys()
    for kk in parentCP2[k].dist.keys():
       if not kk== merger[kk] :
          print '     ', kk ,'merges with ', merger[kk] 
#    if k in parentCP3.keys():
#        print k, 'is caused by:', parentCP3[k].dist.keys()
print 'save data'
f = open(dir+EXPID+'/output/cp/Parentstest.save', 'wb')
cPickle.dump(parentCP2, f, protocol=cPickle.HIGHEST_PROTOCOL)
f.close()
f = open(dir+EXPID+'/output/cp/Kidstest.save', 'wb')
cPickle.dump(kidsCP2, f, protocol=cPickle.HIGHEST_PROTOCOL)
f.close()
#

print 'to loop trough all keys:', sumif2
print 'to loop trough xylist:', sumif3
print 'to loop trough dicton:', sumif1
