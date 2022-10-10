# Why are cell populations maintained via multiple compartments?
# Figure 11
# distribution of G from Gillespie realisations
import numpy as np
import pylab as pl
from random import choice,random
import copy

N = 10
nreal = 10000
print(nreal,'realisations',N)

class Cell(object):
   ''' Cell(g) creates a cell in generation g, compartment c.
   attribute g = generation number '''
   number = 0
   def __init__(self,thisg,thisc):
       self.g = thisg
       self.c = thisc
   def advance(self):
        self.g += 1
   def migrate(self):
        self.c += 1

def onerealization(C,pb,pd):
    ''' Gillespie algorithm, birth-death-differentiation process, 
    C compartments with identical pb and pd'''
    t,livecelllist,deadcelllist,exitcelllist = 0,[Cell(0,1)],[],[]    

    while len(livecelllist)>0:
       thiscell = choice(livecelllist)  # choose a live cell at random
       urv = random()
       if urv < pb:  # division
          thiscell.advance()
          newcell = copy.deepcopy(thiscell)  # copies cell and its attributes
          livecelllist.append(newcell)
       elif urv < 1-pd:  # migration to next compartment
          thiscell.migrate()
          if thiscell.c > C: # exit last compartment
             exitcelllist.append(thiscell)
             livecelllist.remove(thiscell)             
       else: # death
          deadcelllist.append(thiscell)
          livecelllist.remove(thiscell)
    return exitcelllist

def DW(C):
    ''' find D and W given C '''
    i = nreal
    glist = []
    while i > 0:
       prodcelllist = onerealization(C,pb,pd)
       glist.extend([c.g for c in prodcelllist])
       i -= 1
    D = np.mean(glist)
    W = np.var(glist)

    return D,W

#Dfig,(ax1,ax2) = pl.subplots(nrows=2)

import matplotlib
Dfig,(ax1,ax2) = pl.subplots(ncols=2)
ax1.set_yscale('log')
#fig = pl.figure(figsize=(8,3))
params = {'font.family': 'serif'}
matplotlib.rcParams.update(params)

ccolors = pl.rcParams['axes.prop_cycle'].by_key()['color']
myclist = (3,2) 

pd = 0.05
for ic,C in enumerate(myclist):
   pb = (N**(1./C)-1+pd)/(2*N**(1./C)-1)
   print('C = '+str(C))
   i = nreal
   rlist,glist = [],[]
   while i > 0:
      prodcelllist = onerealization(C,pb,pd)
      glist.extend([c.g for c in prodcelllist])
      rlist.append(len(prodcelllist))
      i -= 1
   ax1.hist(rlist,density=True,bins=range(10*N),histtype='step',color=ccolors[ic])
   ax2.hist(glist,density=True,bins=range(N*N),histtype='step',
           linestyle='dotted',color=ccolors[ic])

ax2.set_xlim([0,N*N/2])
ax2.set_ylabel('$\mathbf{P}(\mathbf{G}=n)$',fontsize=14)
ax2.set_xlabel('$n$',fontsize=14)
ax2.set_yticks([0.01,0.02,0.03])
ax1.set_ylabel('$\mathbf{P}(\mathbf{R}=n)$',fontsize=14)
ax1.legend(fontsize=12)
pl.tight_layout()
#pl.savefig('Gdist'+str(N)+'_'+str(int(nreal/1000))+'.jpg',dpi=200)
pl.show()
