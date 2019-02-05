import numpy as np
import numpy.ma as ma
import os   # for environment
import math
import matplotlib
import operator
import matplotlib.pyplot as plt; plt.rcdefaults()
import matplotlib.pyplot as plt
import datetime
from classes import TRACER,CP_map,  CPstart
from operator import itemgetter
from netCDF4 import Dataset
from six.moves import cPickle   # to save class files

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

starttime = datetime.datetime.now()
EXPID ='test1plus4K_fin4'


# SETTINGS
dir = '/nbi/ac/conv1/henneb/results/coldpool/'
odir = '/nbi/home/henneb/Dropbox/CPplots/'  #$HOME

f = open(dir+EXPID+'/output/cp/Parentstest.save', 'rb')
CPparents = cPickle.load(f)
f.close()

f = open(dir+EXPID+'/output/cp/Kidstest.save', 'rb')
CPkids = cPickle.load(f)
f.close()

f = open(dir+EXPID+'/output/cp/CPstart.save', 'rb')
CPstart = cPickle.load(f)
f.close()
start = 1000
end=0
for k in CPstart.keys():
  start= min(start,CPstart[k].time)
  end = max(end,CPstart[k].time)
opentime = datetime.datetime.now()
print 'time to load files:', opentime-starttime

##################################################################
# make a histogram
####
data = {}
for k in CPkids.keys():
    l =  len(CPkids[k].dist.keys())
    if l > 4:
     print k , 'was caused by', CPkids[k].dist.keys()
    if not l in data.keys():
      data[l] = 0
    data[l] += 1
    #print k , 'caused by:', CPkids[k].dist.keys()
for l in data.keys():
   if l > 4:
     print l, data[l]
width = .8

plt.bar(data.keys(), data.values(), width, color="blue")
plt.xticks([a + 0.5 for a in data], data.keys())
plt.show()
plt.savefig(odir+EXPID+'_triggering_hist.png')
#
#####################################################################
data2={}
data2[0]={}
data2[1]={}
data2[2]={}
data2[3]={}
data2[4]={}
datacdf = data2
datacdfrel= data2
for i in range(0,5):
  for j in range(start,end+1):
    data2[i][j] =0

for k in CPparents.keys():
  l =  len(CPparents[k].dist.keys())
  data2[l][CPstart[k].time] +=1 


for i in range(0,5):
  datacdf[i][start]=data2[i][start]
  for j in range(start+1,end+1):
    datacdf[i][j] = datacdf[i][j-1]+data2[i][j]
    


fig,ax = plt.subplots(2,2)
#ax.set_xlim(0, ngp)
#ax.set_ylim(0, ngp)
#plt.xticks(range(0,320,10))
#plt.yticks(range(0,320,10))
print datacdf[0].values()
print datacdf[0][end]
for i in range(0,5):
  ax[0,0].plot(datacdf[i].keys(),datacdf[i].values())
  ax[0,1].plot(datacdf[i].keys(),np.true_divide(list(datacdf[i].values()),datacdf[i][end]))
  

ax[0,1].legend(['0','1','2','3','4'])

plt.show()
