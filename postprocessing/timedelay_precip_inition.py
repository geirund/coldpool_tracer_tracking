import numpy as np
import numpy.ma as ma
import os   # for environment 
import math
import matplotlib
import random
import operator
import matplotlib.pyplot as plt
import datetime
from classes import RAINCLOUD
from operator import itemgetter
from netCDF4 import Dataset
from six.moves import cPickle   # to save class files
EXPID ='test1plus4K' #'lind_p2K' #'test1plus4K' #'lind_p2K' #'test1plus4K'# 'lindp2K'#
DIR = '/nbi/ac/conv1/henneb/results/'
child ='r_int1_10' #'r_int1_10'
v ='w5_4' # 'twp05_20' #'twp05_20' # _50_theta1' #_4_theta0'] #,'rwp1'] # ['rwpa2000m05','rwpb2000m05','twp05', 'rwp05','gwp05','cwp05']#,'gwp05']#,'gwp10']#,'cwp']

data_in_precip_mask    = Dataset(DIR+'/'+EXPID+'/'+child+'/output/irt_tracks_mask.nc', 'r')

# extracting the dimensions of the datafile
dim_x         = data_in_precip_mask.dimensions['x'].size
dim_y         = data_in_precip_mask.dimensions['y'].size
dim_tm        = data_in_precip_mask.dimensions['time'].size

datakid = {}
dataparent = {}

##############################################################
# Variable of last state of rain event (eg sfc precipitation) classified child
# get rain water values for CELLS  from raintracking, rain ID as key, 
# inner key is timestep
# stored info per timestep: age, size, intensity
# stored info valid for all timesteps: start time of track
#############################################################
f = open(DIR+EXPID+'/'+child+'/output/irt_tracks_output_pure_sort.txt', 'r')
lines = f.readlines()
column =lines[0].split()
print column
start_c = int(column[1])   #get the first timestep when tracking starts
print 'precip starts at timestep', start_c
column = lines[-1].split()
end_c = int(column[1])
print 'last timestep in precip tracking', end_c
for line in lines:
    columns = line.split()
    tist = (int(columns[1]))    # timestep
    rainID = (int(columns[0]))
    precip = (float(columns[5]))
    size = (float(columns[4]))
    COGx =(float(columns[8]))
    COGy = (float(columns[9]))
    if not rainID in datakid.keys():
      datakid[rainID]= RAINCLOUD(rainID)
    datakid[rainID].add_info(tist,size,precip)
f.close()

###############################################################
# get rain water values for CELLS from raintracking, rain ID as key,
# inner key is timestep
# stored info per timestep: age, size, intensity
# stored info valid for all timesteps: start time of track
#############################################################
f = open(DIR+EXPID+'/'+v+'/output/irt_tracks_output_pure_sort.txt', 'r')
lines = f.readlines()
column = lines[0].split()
start_p = int(column[1])
print v ,'tracking starts at', start_p
column = lines[-1].split()
end_p = int(column[1])
#print v, 'tracking ends at', end_p

for line in lines:
    columns = line.split()
    tist = (int(columns[1]))    # timestep
    rainID = (int(columns[0]))
    precip = (float(columns[5]))
    size = (float(columns[4]))
    COGx =(float(columns[8]))
    COGy = (float(columns[9]))
    if not rainID in dataparent.keys():
      dataparent[rainID]= RAINCLOUD(rainID)
    dataparent[rainID].add_info(tist,size,precip)

###############################################################
# combine precip and graupel
# read rain track mask as nc
# note: tracks for different variables (rain, graupel) start at differnt timesteps
# loop trough common tiesteps of both masks
# finds all precip (kid) ID and check for beeing the first timestep of the track
# loops trough IDs of starting precip at current timestep
# makes 0-1 mask where precip object is and multiplies with parent mask.
# all IDs in parent mask are then overlaping with this ID and added as list to precipitation class
#############################################################

data_in_rain_mask    = Dataset(DIR+'/'+EXPID+'/'+v+'/output/irt_tracks_mask.nc', 'r')
dim_var        = data_in_rain_mask.dimensions['time'].size

deltaCP = start_c - start_p # how much is child tracking started behind parent? most likely positiv.
mem = {} # remember which kids- sfc precip IDs were used already to identify first timestep
print 'startc',start_c, 'startp',start_p, deltaCP
for ti in range(max([0,deltaCP]),min([dim_var-deltaCP,dim_tm+deltaCP]),1) : #max(dim_tm,dim_var),1):
  print ti, ' von ' , dim_tm
  if ti-deltaCP < dim_tm:
    data_kmask = np.array(data_in_precip_mask.variables["var1"][ti-deltaCP,:,:])
    current_KID_IDs  = np.unique(data_kmask).astype(int)
  else:
    current_KID_IDs = [0,-1]
  if ti < dim_var:
    data_pmask = np.array(data_in_rain_mask.variables["var1"][ti,:,:])
    current_PARENT_IDs  = np.unique(data_pmask).astype(int)
  else:
    current_PARENT_IDs = [0,-1]
### overlap rain and surface precip
  for IDkid in [iid for iid in current_KID_IDs if iid not in[-1,0]]: # (IDi in current_QR_IDs if IDi not in rain.keys()): #precip ID
    #print 'precip ID',IDkid, mem.keys()
    if IDkid not in mem and IDkid not in [0,-1]: # check only first timestep of precip
     mem[IDkid] = IDkid
     mask = data_kmask == IDkid   # mask with only current kid object
     #prec[IDkid]= RAINCLOUD(IDkid)
     data_maskmask = (data_pmask * mask).astype(int)
     rwpID = np.unique(data_maskmask).astype(int)
     #print 'found parents', rwpID
     no_i = 0
     #[[x,l.count(x)] for x in set(l)]
     main_parent=0
     for i in [j for j in rwpID if j not in [-1,0]]:
         #print 'in', i ,'found',(data_maskmask == i).sum()
         no_i_new=(data_maskmask == i).sum() # number ofoverlapping gp
         if no_i_new > no_i:
             main_parent=i
             no_i=no_i_new
     #    print main_parent
     #if len(rwpID) > 0:
     datakid[IDkid].add_parent(v,[i for i in rwpID if i not in [0,-1]],main_parent)

     print('at', ti+1, 'for precip', IDkid, ' found', rwpID)
     data_maskmask = 0
     mask = 0
     main_parent = 0


#################
# get timedelay between time when graupel object starts and when precip object starts
# read delay file for timedifference between start of precip track and time when tracers are set
# which is either when track size exceed treshold value or for tracks smaller than treshold when they reach their max size
# if several graupel objects overlap with precip possible decission:
#     1. take largest/ most intense (peak intense or sum)
#     2. take earliest
# and get maximal precip intensity
#

# read delay file
f = open('/nbi/ac/conv1/henneb/results/coldpool/test1plus4K/input/cp/mergingCPs.txt', 'r')
lines = f.readlines()
for line in lines:
    columns = line.split()
    rID = (int(columns[0]))
    trackstart = (int(columns[3]))
    tracerstart = (int(columns[2]))
    datakid[rID].trackdelay = tracerstart-trackstart
    mergeID = (int(columns[1]))
    #if mergeID == 0:
     # del datakid[rID]
for k in [kk for kk in datakid.keys() if not kk in [-1,0]]:
  if not datakid[k].start == min(datakid[k].ts.keys()):
           print 'wrong start', datakid[k].start , min(datakid[k].ts.keys())
  #print datakid[k].parent.keys()
## variante parent with highest value
#  if v in datakid[k].parent.keys():
#    #print 'length  of kids:', len(datakid[k].parent[v])
#    lenkid =len(datakid[k].parent[v])
#    if lenkid == 0:
#      print 'no associated cloud',v,'found to surface precipitation with ID', k
#    elif lenkid > 1:
#      #print 'going into routine of more than 1 kid'
#      int_temp = 0
#      # loop trough all parent IDs connected with this cell
#      for parID in datakid[k].parent[v]:
#         # chose the strongest parent
#         # chose the parent with the highes maximal(over time) intensity
#         maxi = max(dataparent[parID].inten.values())
#         if maxi > int_temp:
#           parentID = parID
#           int_temp = maxi
#         # timedifference between start of graupel and start of sfc precip
#         datakid[k].deltat(datakid[k].start - dataparent[parentID].start)
#         print 'precipitation',k,'starting at', datakid[k].start, 'with graupel',parentID,'istarting at',dataparent[parentID].start
#####################################################
#          # time between half way graupel and start of sfc precip
#          # half way graupel = graupel at half its max value
#         tp = min(dataparent[parentID].inten.keys())
#         while dataparent[parentID].inten[tp] < dataparent[parentID].inten[datakid[k].start]/2.:
#              tp+=1
#         #datakid[k].deltat(max(datakid[k].inten, key=datakid[k].inten.get) - datakid[k].start ) #tp)
#         datakid[k].deltat(datakid[k].start  -tp)
#    else:
#        parentID = datakid[k].parent[v][0]
##       # timedifference between start of graupel and start of sfc precip
##        datakid[k].deltat(datakid[k].start - dataparent[parentID].start)
######################################################3
#        # time between half way graupel and start of sfc precip
#        # half way graupel = graupel at half its max value
#        tp = min(dataparent[parentID].inten.keys())
#        while dataparent[parentID].inten[tp] < dataparent[parentID].inten[datakid[k].start]/2.:
#             tp+=1
#        datakid[k].deltat(datakid[k].start  -tp)
##       datakid[k].deltat( max(datakid[k].inten, key=datakid[k].inten.get)- datakid[k].start ) #tp)
########################################################33
# Variante with most overlap
  parentID  = datakid[k].mainparent
  #print 'k',k, 'parent', parentID
  if not parentID ==0:  # default if not parent was found
      ttemp = datakid[k].start
      print 'kidtimes', datakid[k].inten.keys()
      print 'kid', k, 'parent', parentID, ttemp
      print dataparent[parentID].inten.keys()
      tempinten = dataparent[parentID].inten[ttemp]
      while dataparent[parentID].inten[ttemp] <= tempinten:
          tstart = ttemp
          if ttemp-1 in dataparent[parentID].inten.keys():
              ttemp = ttemp-1
          else:
              tempinten = 0

      datakid[k].deltat(datakid[k].start -tstart)   #dataparent[parentID].start)
#      datakid[k].deltat(datakid[k].start -dataparent[parentID].start)

#######################################
# make scatter plot
#####################################
avgdelay={}
count = {}
avgtime={}
fig, ax = plt.subplots(1)  # make some smart devision accoridng to number of varaibles
for k in datakid.keys():
    td = datakid[k].delaycloud
    if not td in avgdelay.keys():
        avgdelay[td] = max(datakid[k].inten.values())
        avgtime[td]=(float(datakid[k].start)-float(start_c))/float((end_c-start_c))
        count[td] = 1
    else:
        avgdelay[td]+= max(datakid[k].inten.values())
        count[td] +=1
        avgtime[td]+=(float(datakid[k].start)-float(start_c))/float((end_c-start_c))
for k in avgdelay.keys():
    avgdelay[k] = float(avgdelay[k])/float(count[k])
    avgtime[k] = float(avgtime[k])/float(count[k])
plotdatax = {}
plotdatay = {}
plotdatac = {}
for k in datakid.keys():
  plotdatay[k] = datakid[k].delaycloud
  plotdatax[k] = max(datakid[k].inten.values())
  plotdatac[k] = (float(datakid[k].start)-float(start_c))/float((end_c-start_c))

#for k in plotdatax.keys():
ax.scatter(plotdatax.values(),plotdatay.values(),s=30, cmap='jet',c=plotdatac.values(), facecolors=plotdatac.values(),edgecolors="None",alpha=0.75) #"None")
  #ax.scatter(plotdatax[k],plotdatay[k],s=60, cmap='hsv',facecolor=plotdatac[k],edgecolor="None")
#ax.scatter(max(datakid[k].inten.values()),datakid[k].delay,facecolors=(float(datakid[k].start)-float(start_p))/float((end_p-start_p))) #,edgecolors="None")
#  ax[j].set_aspect('equal', 'box')
ax.set_xlabel('max rain intensity/ mm'+'$\mathregular{h^{-1}}$')
ax.set_ylabel('time delay/ (5 min)')
#ax.title.set_text(v)
ax.set_ylim(-0.5,13.5)
im=ax.scatter(avgdelay.values(),avgdelay.keys(),s=130,cmap='jet',c=avgtime.values(),edgecolors="None")
fig.colorbar(im, ax=ax)
print avgdelay.values(), avgdelay.keys()
for i in datakid.keys(): #[143,148,163,165,166,193]:
  print i, datakid[i].delay, 'overlap', datakid[i].mainparent,datakid[i].start-136

f = open(DIR+EXPID+'/r_int1/output/timedelay.save', 'wb')
cPickle.dump(datakid, f, protocol=cPickle.HIGHEST_PROTOCOL)
f.close()

plt.show()
fig.savefig('/nbi/home/henneb/plots/'+EXPID+'/timedelay/maxprec_startcloudtostartprec'+v+'.png')

