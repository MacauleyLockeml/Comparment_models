# Why are cell populations maintained via multiple compartments?
# Figure 6
# Q_k(C) vs k, multiply by k^(3/2)
# one and multiple compartments
# direct histogram
import numpy as np
import pylab as pl
from random import choice,random
import copy
from collections import defaultdict,Counter

import matplotlib
fig = pl.figure(figsize=(8,3))
params = {'font.family': 'serif'}
matplotlib.rcParams.update(params)

N = 25
nreal = 300000
mypd = 0.0
filename = 'Qkhist'+str(N)+'_'+str(nreal/10000)+'.png'

print([float("%0.2f"%i) for i in (N,nreal,mypd)])

def qasympt(k):
    '''asymptotic formula C=1'''
    delta = np.sqrt(1-4*pb*pd)
    if k==0:
        return 2*pd/(1+delta)
    elif k==1:
        return pe/delta
    else:
        fac = (4*pb*pe)/(1-4*pb*pd)
        return pe/delta*fac**(k-1)/((k-1)*np.sqrt(np.pi*(k-1)))

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
    return len(exitcelllist)

pl.subplot(111)

mycol = {1:'r',2:'b',3:'g',4:'y',10:'c'}

for C in (1,2,10):
    pd = mypd/C
    pb = (N**(1./C)-1+pd)/(2*N**(1./C)-1)
    pe = 1-pd-pb
    delta = np.sqrt(1-4*pb*pd)
    print([float("%0.4f"%i) for i in (C,pb,pd,pe,N**(1./C))])
    i = nreal
    qlist = []
    while i>0:
        qlist.append(onerealization(C,pb,pd))
        i -= 1
        if i==nreal/2:
            print('halfway')
    qhist = Counter(qlist)
    n = 64
    klist,lhist = [n],[qhist[n]*1./nreal]
    while n < max(qhist.keys()):
        w = int(np.sqrt(n*0.4))
        thisk = n
        lhist.append(thisk**1.5*sum([qhist[j] for j in range(n-w,n+w)])*1./(2*w*nreal))
        klist.append(thisk)
        n = int(n*1.25)
    pl.loglog(klist,lhist,'o',color=mycol[C],label='C='+str(C))
    a2 = 1 - 0.25/N**(1+1./C)
#    pl.loglog([4*N**(1+1./C),4*N**(1+1./C)],[0.01,1],color=mycol[C],linestyle=':')

    if C==1:
        q0 = 2*pd/(1+delta)
        q1 = pe/delta
        fac = 2*pb*pe/delta**2
        cat = [q0,q1]
        n=1
        while n < 20000:
            cat.append(fac*cat[-1]*(2*n-1.)/(n+1))
            n += 1
        pl.loglog([q*i**1.5 for i,q in enumerate(cat)],color=mycol[C])
    if C==2:
        Qk = [0,pe*pe,pb*(pe**3+pe**4)]
        k=1
        while k < 3500:
            f1 = (1+1./k)*(1+2./k)*(1-4*pe)
            f2 = (2+1./k)*(2+2./k)*pb*pe*(1-4*pe-4*pe*pe)
            f3 = (64-4./k**2)*pb**2*pe**4
            thisQ = (f2*Qk[-1] + f3*Qk[-2])/f1
            Qk.append(thisQ)
            k += 1
        pl.loglog([q*i**1.5 for i,q in enumerate(Qk)],color=mycol[C])

print([(klist[i],lhist[i]) for i in range(len(klist)) if klist[i]>300])

pl.xlim([60,8000])
pl.xlabel('$k$',fontsize=14)
#pl.ylim([2e-4,3.])
pl.ylim([0.008,2.])
#pl.yticks([1e-6,1e-4,1e-2],fontsize=10)
pl.ylabel('$k^{3/2}Q_k(C)$',fontsize=14)
pl.legend(fontsize=12)
pl.tight_layout()
pl.savefig(filename,dpi=200)
pl.show()


