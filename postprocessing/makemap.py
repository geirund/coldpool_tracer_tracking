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
EXPID ='test1plus4K_fin4' # 'lindp2K_fin' #'test1plus4K_fin2'
DIR = '/nbi/ac/conv1/henneb/results/test1plus4K/' #lindp2K/' #test1plus4K/'
# SETTINGS
dir = '/nbi/ac/conv1/henneb/results/coldpool/'
ngp = 320
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
cmap = plt.get_cmap('jet') #'terrain') #viridis')
def arrowsetting(x1,x2,y1,y2,clr):
  back = 5  # how much shell arrow end before point
  DELTAx = x2-x1
  DELTAy = y2-y1
  alpha = math.atan2(DELTAy,DELTAx)
  print alpha
  dx = back * np.cos(alpha)
  dy = back * np.sin(alpha)
  ax.arrow(x1,y1,DELTAx-dx,DELTAy-dy,head_width=3,head_length=2, zorder=2,color=clr)


def set(i,CPinit,kids,c):
  x  = CPinit[i].x
  y  = CPinit[i].y
  outc[i] =CPinit[i].time # (CPinit[i].time - starttime) *2
  #ax.plot(x,y,'o', markersize=15,markerfacecolor=cmap(outc[i]*40),markeredgecolor="None",zorder=1)
  #ax.plot(x,y,'o', markersize=15,markerfacecolor=cmap(outc[i]),markeredgecolor="None",zorder=1)
   
  outx[i] = CPinit[i].x
  outy[i] = CPinit[i].y
  print i, c
  if i in kids.keys():
    c += 1
    nkids= len(kids[i].dist.keys())
#    print nkids
    if nkids > 0 :
#      print kids[i].dist.keys()
      for j in kids[i].dist.keys():
        if j in outc.keys():
          outc[j] = max(outc[i] +1, outc[j])
        else:
          outc[j] = outc[i] +1
        dx = CPinit[j].x - CPinit[i].x 
        dy = CPinit[j].y - CPinit[i].y

        if not kids[i].dist[j] > 100:
          #ax.plot([CPinit[i].x,CPinit[j].x], [CPinit[i].y, CPinit[j].y], 'k-', lw=2)
          arrowsetting(CPinit[i].x,CPinit[j].x,CPinit[i].y,CPinit[j].y,'black') 
        elif dx < -160 and dy > -160 and dy < 160:
          arrowsetting(CPinit[i].x-ngp,CPinit[j].x,CPinit[i].y,CPinit[j].y,'black')
          arrowsetting(CPinit[i].x,CPinit[j].x+ngp,CPinit[i].y,CPinit[j].y,'black')
        elif dx > 160 and dy > -160 and dy < 160:
          arrowsetting(CPinit[i].x+ngp,CPinit[j].x,CPinit[i].y,CPinit[j].y,'black')
          arrowsetting(CPinit[i].x,CPinit[j].x-ngp,CPinit[i].y,CPinit[j].y,'black')
        elif dy < -160 and dx > -160 and dx < 160:
          arrowsetting(CPinit[i].x,CPinit[j].x,CPinit[i].y-ngp,CPinit[j].y,'black')
          arrowsetting(CPinit[i].x,CPinit[j].x,CPinit[i].y,CPinit[j].y+ngp,'black')
        elif dy > 160 and dx > -160 and dx<160:
          arrowsetting(CPinit[i].x,CPinit[j].x,CPinit[i].y+ngp,CPinit[j].y,'black')
          arrowsetting(CPinit[i].x,CPinit[j].x,CPinit[i].y,CPinit[j].y-ngp,'black')

        set(j,CPinit,kids,c)

f = open(dir+EXPID+'/output/cp/CPstart.save', 'rb')
CPinit = cPickle.load(f)
f.close()

starttime = 300
for k in CPinit.keys():
 starttime = min(CPinit[k].time, starttime)


f = open(dir+EXPID+'/output/cp/Parentstest.save', 'rb')
parents = cPickle.load(f)
f.close()

f = open(dir+EXPID+'/output/cp/Kidstest.save', 'rb')
kids = cPickle.load(f)
f.close()

outx = {}
outy= {}
outc = {}

fig, ax = plt.subplots(1)
ax.set_xlim(0, ngp)
ax.set_ylim(0, ngp)
plt.xticks(range(0,320,10))
plt.yticks(range(0,320,10))
#ax.xaxis.set_major_locator(MaxNLocator(4))
# get merger
f = open(dir+EXPID+'/input/cp/mergingCPs.txt', 'rb')
lines = f.readlines()
#for line in lines:
#  column=line.split()
#  raincell= (int(column[0])) 
#  rainmerger= (int(column[1]))
#  if rainmerger != 0:
#    if rainmerger != raincell:
#      if not ((CPinit[rainmerger].x-CPinit[raincell].x)**2+ (CPinit[rainmerger].y- CPinit[raincell].y)**2)**(0.5) > 160:
#       #ax.plot([CPinit[rainmerger].x,CPinit[raincell].x], [CPinit[rainmerger].y, CPinit[raincell].y], 'k-', lw=5,c='red')
#       ax.arrow(CPinit[rainmerger].x,
#                CPinit[rainmerger].y, 
#                CPinit[raincell].x-CPinit[rainmerger].x,
#                CPinit[raincell].y-CPinit[rainmerger].y, color='red')
        
#print kids.keys
c = 0
for i in CPinit.keys():
  nparents = len(parents[i].dist.keys())

  print nparents
  if nparents == 0:
    outc[i] = 0
    print 'at',i
    set(i,CPinit,kids,c)

f = open(dir+EXPID+'/input/cp/mergingCPs.txt', 'rb')
lines = f.readlines()

for line in lines:
  column=line.split()
  raincell= (int(column[0]))
  rainmerger= (int(column[1]))
  if raincell in CPinit.keys() and rainmerger in CPinit.keys():
   if rainmerger != 0:
    if rainmerger != raincell:
      if not ((CPinit[rainmerger].x-CPinit[raincell].x)**2+ (CPinit[rainmerger].y- CPinit[raincell].y)**2)**(0.5) > 160:
       arrowsetting(CPinit[rainmerger].x,CPinit[raincell].x,CPinit[rainmerger].y,CPinit[raincell].y,'red')
        
im = ax.scatter(outx.values(), outy.values(), c=outc.values(),s=90,edgecolor='none',cmap=cmap)
fig.colorbar(im, ax=ax)
#cbar = fig.colorbar(im, ticks=[-1, 0, 1])
plt.show()
#fig.savefig(odir+'/plots/map_rain_tracer'+str(data.ts-1)+'.png')
   
