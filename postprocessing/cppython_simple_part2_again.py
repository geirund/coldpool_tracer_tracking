import numpy as np
import os   # for environment 
import datetime
from classes import CP_map, CPlife
import matplotlib.pyplot as plt
from six.moves import cPickle 


EXPID = 'lindp2K_fin'#
atime = datetime.datetime.now()
odir = os.environ.get('results')+'/coldpool/'
f = open(odir+EXPID+'/output/cp/TracerMap.save', 'rb')
CoPo = cPickle.load(f)
f.close()

b = {}
for (t,x,y) in CoPo.keys() :
    #print x,y,t     
    cCP = CoPo[t,x,y].CPs
    #print 'cps', list(cCP.values())
#    # summerize collisions:
#    if zip(*zip(sorted(cCP.values()))[0]) not in CPL.keys():
#        print list(cCP.values()), CPL.keys()
#        #if this CP combination not yet exists, make it
#        CPL[zip(*zip(sorted(list(cCP.values()))))[0]]= COL()
#    else:
#        print x,y,t, 'list is added'
#    # add location of tracers 
#    CPL[zip(*zip(sorted(list(cCP.values()))))[0]].add(t,x,y,CoPo[t,x,y].nTrtot)

    # summerize CP lifecycle
    for cp in list(cCP.values()):
      if cp not in b.keys():
        b[cp] = CPlife(cp,t)
      if len(list(cCP.values())) ==1 :
        b[cp].add(CoPo[t,x,y].nTrCP[cp],CoPo[t,x,y].nTrCP[cp],1,t,CoPo[t,x,y].age[cp],x,y)
      else:
        b[cp].add(CoPo[t,x,y].nTrCP[cp],0,0,t,CoPo[t,x,y].age[cp],x,y)
        for cps in list(cCP.values()):
         if not cps == cp:
          b[cp].add_others(t,cps,CoPo[t,x,y].nTrCP[cp],CoPo[t,x,y].nTrCP[cps] )

ctime = datetime.datetime.now()
print 'took ',(ctime-atime), 'for loop '
f = open(odir+EXPID+'/output/cp/CPlife.save', 'wb')
cPickle.dump(b, f, protocol=cPickle.HIGHEST_PROTOCOL)
f.close()

#
