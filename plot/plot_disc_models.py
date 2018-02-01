# Written 19/9/2017 by dh4gan
# Reads in Ken's pre-run disc models and plots them
# Optional - sample where fragments can be found

import numpy as np
import matplotlib.pyplot as plt

filename = 'migrate1_rout100.dat'

nr = 200
AU = 1.496e13
Q=2.0
aspectratio = 0.1

step = 1

t, Lx = [], []
r, sigma, nu, Tc = [], [], [], []
tau, alpha, Mstar = [], [], []

with open(filename,'r') as inputfile:
    for j in range(step+1):
        print j
        fullline = inputfile.readline()
        line = fullline.split()        
        
        t.append(line[0])
        Lx.append(line[1])
        
        for i in range(nr):
            fullline = inputfile.readline()
            line = [float(i) for i in fullline.split()]
            print line
            
            if(j==step):
                r.append(line[0])
                sigma.append(line[1])
                nu.append(line[2])
            
            fullline = inputfile.readline()
            line = [float(i) for i in fullline.split()]
            print line
            
            if(j==step):
                Tc.append(line[0])
                tau.append(line[1])
                alpha.append(line[2])
            
            line = inputfile.readline()
            
            if(j==step):
                Mstar.append(float(line))            
            
r = np.array(r)    
r = r/AU

alpha = np.array(alpha)

# Find where alpha exceeds 0.06 (fragmentation boundary, roughly same as Jeans mass calculation in this setup)


ifrag = np.where(alpha[:]>0.06)[0][0]

print ifrag, r[ifrag]

# Add fragments according to pop synthesis description (read from a model run)


rfrag = np.array([19.2,28.1,51.85,79.69])
#mfrag = np.array([50.85,77.7,87.8,107.0])
mfrag = np.array([100.0,100.0,100.0,100.0])

sigcolor = 'red'
alphacolor = 'blue'


fig1 = plt.figure(figsize=(10,8))
ax1 = fig1.add_subplot(111)
ax1.set_xscale('log')
ax1.set_xlim(0.5,200.0)
ax1.set_ylim(1.0e1,1.0e5)

ax1.set_xlabel('r (AU)', fontsize=25)
ax1.set_yscale('log')
ax1.set_ylabel(r'$\Sigma$ (g cm$^{-2}$)',fontsize=25,color=sigcolor)
ax1.tick_params(labelsize=20)

lfrag = ax1.scatter(rfrag,mfrag,s=40,color='black',marker='o', label='Fragments')

ax2 = ax1.twinx()
ax2.set_xlim(0.5,200.0)
ax2.set_yscale('log')
ax2.set_ylabel(r'$\alpha$',fontsize=25,color=alphacolor)
ax2.tick_params(labelsize=20)

lsig = ax1.plot(r,sigma,linewidth=3,color=sigcolor,label=r'$\Sigma$')
lalpha = ax2.plot(r,alpha,linewidth=3, linestyle='dashed',color=alphacolor,label=r'$\alpha$')

ax1.legend(loc='upper left',fontsize=22)
ax2.legend(loc='lower left',fontsize=22)

plt.show()

fig1.savefig('discprofile.png')
