#
# Written 28/3/17 by dh4gan
# Plots the orbital architecture of output systems
# generated from grapus v3.0
# Selects a random sample from system
#

import matplotlib.pyplot as plt
import matplotlib.patches as patches
import filefinder as ff
import numpy as np
from io_grapus import nfinalcol,finalcoldict,finallabeldict, logcoldict,brown,red,blue

print ''
print '\t \t \t ------------------------'
print '\t \t \t Tidal Downsizing System Portraits (N body data)'
print '\t \t \t ------------------------'
print ''


plt.rcParams['font.size'] = 18
plt.rcParams['patch.linewidth'] = 2
plt.rc('xtick', labelsize=14)
plt.rc('ytick', labelsize=14)

nsystems = 10 # Number of systems to plot
ntrymax = 100000 # Number of samplings to try
systemspacing = 1000
patchwidth = 1000

markerscale = 3
minmarkersize = 50
amin = 0.1
amax=1.0e5


# Use filefinder to find files

finalfile = ff.find_sorted_local_input_files('*.final')

print 'Reading ',finalfile

finaldata = np.genfromtxt(finalfile)
nfinalrow = finaldata.size/nfinalcol
finaldata.reshape(nfinalrow,nfinalcol)

print 'There are initially ',nfinalrow, 'rows'

# Delete rows where semi major axis equal to the boundary value rin

istarcol = finalcoldict['istar']
acol = finalcoldict['a']
ecol = finalcoldict['e']
inccol = finalcoldict['i']
mcol = finalcoldict['mass']
mcorecol = finalcoldict['mcore']

finaldata = finaldata.compress(finaldata[:,acol]>0.11, axis=0)
# Delete low mass bodies without cores

blob = (finaldata[:,mcol]<1.0e-2) & (finaldata[:,mcorecol]==0.0)

finaldata = finaldata[np.logical_not(blob)]

nfinalrow = finaldata.shape[0]

print 'After deletion at inner system boundary, there are ',nfinalrow, ' rows'

print 'Data read'

# Separate data into ejected and non-ejected bodies

ejectadata = finaldata[finaldata[:,ecol]>=1.0]
nejected = ejectadata.shape[0]

finaldata = finaldata[finaldata[:,ecol]<1.0]
nbound = finaldata.shape[0]

isystem = 0
print 'There are ',nejected, ' ejected bodies'
print 'There are ',nbound, ' bodies still bound to their parent stars'

imax = int(np.amax(finaldata[:,istarcol]))

print 'Have ', imax, ' systems to choose from for portraits'

fig1 = plt.figure()
ax1 = fig1.add_subplot(111)

names = ('Gas Giant, \nIcy Core', 'Gas Giant, \nRocky Core','Gas Giant, \nNo Core', 'Terrestrial Planet', 'Brown Dwarfs', 'Planetesimals')

for i in range(len(names)):
    print i, names[i]

types = np.zeros((len(names)))

# (imelt ivap idiss igrown iself ijeans itidal)

# Objects with icy cores
# (0 0 0 1 1 1 *)  in binary = (56,120)

icy = np.array([56,120])

# Objects with rocky cores
# (1 * * 1 1 1 *) in binary = (57,59,61,63,121,122,125,127)

rocky = np.array([57,59,61,63,121,122,125,127])

# Brown Dwarfs/ objects without any core
# (1 1 * * * 0 *) in binary =  (3,7,11,15,19,23,27,31,68,71,75,79,83,87,91)
nocore = np.array([3,7,11,15,19,23,27,31,68,71,75,79,83,87,91])

# Objects which produce planetesimal belts
#(* 0 0 * 0 0 1) in binary = (72,73)
 
planetesimal = np.array([72,73])

alltypes = np.concatenate((icy,rocky,nocore,planetesimal),axis=0)
counter = 0

selected = []

ntries = 0
while isystem < nsystems:

    iselect = np.random.randint(0,high=imax)
    
    ntries = ntries+1
    if(ntries>ntrymax): break
    # Find first instance of iselect

    output = finaldata[finaldata[:,istarcol]==iselect,:]

    #print output[:,mcorecol]*0.003146/output[:,mcol]
    nplanets = output.shape[0]
    
    #if(nplanets<2):continue
    #if(np.amax(output[:,mcorecol]*0.003146/output[:,mcol])<0.5 or (np.amin(output[:,acol])<0.1 or np.amax(output[:,acol])>amax)): continue
    #if(np.amax(output[:,mcol])/np.amin(output[:,mcol])<50.0): continue

    # Only plot systems with more than 3 bodies at a<10 AU
    if(nplanets<2 or np.amin(output[:,acol])>10.0): continue
    
    
    # Check system hasn't been picked before
    if iselect in selected: continue
    
    print 'Selected system ', iselect
    selected.append(iselect)
    isystem = isystem+1

    yposition = np.zeros(nplanets)
    yposition[:]=isystem*systemspacing

    # Marker sizes
    sizes = markerscale*output[:,mcol]

    sizes[sizes<minmarkersize] = minmarkersize
    print sizes
    colors = []
   
        
    print "Classifying Fragments"

    for i in range(nplanets):

        # Convert data into binary code
        binary = 0
        for j in range(1,7):
            binary += int(output[i,j])*2**(j-1)    
    
        # Check for icy core
        #print binary,np.sum(np.argwhere(icy==binary))>0,np.sum(np.argwhere(rocky==binary))>0,np.sum(np.argwhere(nocore==binary))>0,np.sum(np.argwhere(planetesimal==binary))>0 
    
        # Icy Core
        if(np.sum(np.argwhere(icy==binary)>0)): colors.append(red)
    
        # Rocky Core
        if(np.sum(np.argwhere(rocky==binary)>0)):
        
            # Is this gas giant or terrestrial planet? determine by mcore/membryo
            mratio = output[i,mcorecol]*0.003146/output[i,mcol]
         
            if(mratio < 0.5):
                colors.append(red)
            else:
                colors.append(blue)
    
        # No Core /BD
        if(np.sum(np.argwhere(nocore==binary)>0)):
        # BD 
            if output[i,mcol] >13.0 and output[i,3]==1:
                colors.append(brown)
            # No core
            else:
                colors.append(blue)
                
    print output[:,mcol],output[:,acol],output[:,ecol]
    ax1.add_patch(patches.Rectangle( (amin,yposition[0]-0.5*patchwidth), amax, patchwidth, color='gray', alpha=0.2))
    ax1.scatter(output[:,acol],yposition,s=sizes, color=colors)
    ax1.errorbar(output[:,acol],yposition,xerr=output[:,acol]*output[:,ecol], linestyle='None', color='black')
    

ax1.set_xlim(1.0e-1,amax)
ax1.set_xscale('log')
ax1.set_ylim(0,yposition[0]+patchwidth)
ax1.set_xlabel('a (AU)',fontsize=22)
ax1.get_yaxis().set_ticks([])
plt.show()
fig1.savefig('portrait.png')

