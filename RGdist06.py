# Why are cell populations maintained via multiple compartments?
# Figure 12
# distributions of R and G from Gillespie realisations
# model where first compartment has asymmetric division
import numpy as np
import pylab as pl
from random import choice,random
import copy

C = 5
pd1list = (0.55,0.90)  # list of death rates in first compartment 
pd2, pd4 = 0.45,0.25  # death rates in subsequent compartments 
pe = 0.3  # values for c=2,3,4,5

nreal = 20000000
print(nreal,'realisations','p_e =',pe)

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
    t,livecelllist,deadcelllist,exitcelllist = 0,[Cell(0,1)],[]    ,[]    

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

def onerealizationa(C,pb,pd,pd1):
    ''' Gillespie algorithm, birth-death-differentiation process, 
    only p_a and p_d \ne 0 in first compartment, 
    then C-1 compartments with identical parameters and p_a=0'''
    t,livecelllist,earlycelllist,exitcelllist = 0,[],[Cell(0,1)],[]    

    sumlen = len(livecelllist)+len(earlycelllist)
    while sumlen > 0:
       urv = random()
       if random() < len(earlycelllist)*1./sumlen:
          thiscell = choice(earlycelllist)  # choose a cell with c = 1
          if urv < pd1:
             earlycelllist.remove(thiscell)
          else:
             thiscell.advance()
             newcell = copy.deepcopy(thiscell)
             newcell.migrate()
             livecelllist.append(newcell) 
       else:
          thiscell = choice(livecelllist)  # choose a cell with c > 1
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
             livecelllist.remove(thiscell)
       sumlen = len(livecelllist)+len(earlycelllist)
    return exitcelllist

def onerealizationb(C,pd1,pd2,pd4):
    ''' Gillespie algorithm, birth-death-differentiation process, 
    only p_a and p_d \ne 0 in first compartment, 
    then two compartments with N_i < 1
    and two more compartments with N_i>1.
    p_b(2) = 1-pe-p_d(2)  and  p_b(4) = 1-pe-p_d(4)'''
    t,latecelllist,midcelllist,earlycelllist,exitcelllist = 0,[],[],[Cell(0,1)],[]    

    sumlen = len(latecelllist)+len(midcelllist)+len(earlycelllist)
    while sumlen > 0:
       urv1, urv2 = random(),random()
       if urv1 < len(earlycelllist)*1./sumlen:
          thiscell = choice(earlycelllist)  # choose a cell with c = 1

          if urv2 < pd1:
             earlycelllist.remove(thiscell)
          else:
             thiscell.advance()
             newcell = copy.deepcopy(thiscell)
             newcell.migrate()
             midcelllist.append(newcell) 

       elif urv1 > 1- len(latecelllist)*1./sumlen:
          thiscell = choice(latecelllist)  # choose a cell with c > 3

          if urv2 < 1-pe-pd4:  # division
             thiscell.advance()
             newcell = copy.deepcopy(thiscell)  # copies cell and its attributes
             latecelllist.append(newcell)
          elif urv2 < 1-pd4:  # migration to next compartment
             thiscell.migrate()
             if thiscell.c > C: # exit last compartment
                exitcelllist.append(thiscell)
                latecelllist.remove(thiscell)             
          else: # death
             latecelllist.remove(thiscell)
       else:
          thiscell = choice(midcelllist)  # choose a cell with c = 2 or 3

          if urv2 < 1-pe-pd2:  # division
             thiscell.advance()
             newcell = copy.deepcopy(thiscell)  # copies cell and its attributes
             midcelllist.append(newcell)
          elif urv2 < 1-pd2:  # migration to next compartment
             thiscell.migrate()
             if thiscell.c > C-2: # move to next type
                latecelllist.append(thiscell)
                midcelllist.remove(thiscell)             
          else: # death
             midcelllist.remove(thiscell)

       sumlen = len(latecelllist)+len(midcelllist)+len(earlycelllist)
    return exitcelllist

import matplotlib
Dfig,(ax1,ax2) = pl.subplots(nrows=2)
ax1.set_yscale('log')

Dfig.set_size_inches(7,3)
#fig = pl.figure(figsize=(8,3))
params = {'font.family': 'serif'}
matplotlib.rcParams.update(params)

ccolors = pl.rcParams['axes.prop_cycle'].by_key()['color']

for ic,pd1 in enumerate(pd1list):
   print('pd(c) values:',pd1,pd2,pd2,pd4,pd4)
   i = nreal
   rlist,glist = [],[]
   while i > 0:
      prodcelllist = onerealizationb(C,pd1,pd2,pd4)
      glist.extend([c.g for c in prodcelllist])
      rlist.append(len(prodcelllist))
      i -= 1

   print(' R',np.mean(rlist),np.var(rlist))
   print(' G',np.mean(glist),np.var(glist))

   ax1.hist(rlist,density=True,bins=range(160),histtype='step',
            color=ccolors[ic])

   if ic==0:
      ax2.hist(glist,density=True,bins=range(int(110)),rwidth=0.5,
               color=ccolors[ic],align='left',label='$p_a(1)=0.45$')
   else:
      ax2.hist(glist,density=True,bins=range(int(110)),rwidth=0.5,
               color=ccolors[ic],align='mid',label='$p_a(1)=0.10$')

ax1.tick_params(axis='x', labelsize=8)
ax2.tick_params(axis='x', labelsize=8)
ax1.tick_params(axis='y', labelsize=8)
ax2.tick_params(axis='y', labelsize=8)

ax2.set_xlim([0,80])
ax2.set_ylabel('$\mathbf{P}(\mathbf{G}=n)$',fontsize=14)
ax2.set_xlabel('$n$',fontsize=12)
ax2.set_yticks([0.01,0.02,0.03])

ax1.set_xlim([0,130])
ax1.set_xlabel('$n$',fontsize=12)
ax1.set_ylabel('$\mathbf{P}(\mathbf{R}=n)$',fontsize=14)
ax1.set_yticks([0.1,0.001,0.00001])

ax2.legend(fontsize=12)
pl.tight_layout()
pl.savefig('RG6dist'+'_'+str(int(pd1list[0]*10))+'_'+str(int(nreal/1000))+'.jpg',dpi=200)
pl.show()
