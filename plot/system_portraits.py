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
from io_grapus import read_finaldata, get_fragment_colour, finalcoldict

print ''
print '\t \t \t ------------------------'
print '\t \t \t Tidal Downsizing System Portraits (N body data)'
print '\t \t \t ------------------------'
print ''




nsystems = 50 # Number of systems to plot
ntrymax = 100000 # Number of samplings to try
systemspacing = 1000
patchwidth = 0.7*systemspacing
figheight = nsystems
figwidth = 0.8*nsystems

plt.rcParams['font.size'] = 18*np.int(nsystems)/10
plt.rcParams['patch.linewidth'] = 2
plt.rc('xtick', labelsize=14*np.int(nsystems)/10)
plt.rc('ytick', labelsize=14*np.int(nsystems)/10)



markerscale = 3
minmarkersize = 50
amin = 0.1
amax=1.0e5

# Use filefinder to find files

finalfile = ff.find_sorted_local_input_files('*.final')

finaldata, ejectadata, nbound,nejected = read_finaldata(finalfile)

isystem = 0

imax = int(np.amax(finaldata[:,finalcoldict['istar']]))

print 'Have ', imax, ' systems to choose from for portraits'



fig1 = plt.figure(figsize=(figheight,figwidth))
ax1 = fig1.add_subplot(111)

counter = 0

iselected = []

allsystems = []
allypositions = []
allsizes = []
allcolors = []

ntries = 0
while isystem < nsystems:

    iselect = np.random.randint(0,high=imax)
    
    ntries = ntries+1
    if(ntries>ntrymax): break
    # Find first instance of iselect

    output = finaldata[finaldata[:,finalcoldict['istar']]==iselect,:]

    nplanets = output.shape[0]
    if nplanets<1: continue
    
    # Only plot systems with more than 3 bodies at a<10 AU
    
    #if(nplanets<2 or np.amin(output[:,finalcoldict['a']])>10.0): continue
    
    
    # Check system hasn't been picked before
    if iselect in iselected: continue
    
    print 'Selected system ', iselect
    iselected.append(iselect)
    allsystems.append(output)
    
    isystem = isystem+1

    yposition = np.zeros(nplanets)
    yposition[:]=isystem*systemspacing
    
    allypositions.append(yposition)

    # Marker sizes
    sizes = markerscale*output[:,finalcoldict['mass']]

    sizes[sizes<minmarkersize] = minmarkersize

    allsizes.append(sizes)

    colors = []
   
   # print "Classifying Fragments"

    for i in range(nplanets):
        colors.append(get_fragment_colour(output,i))
    
    allcolors.append(colors)

print allypositions
for iselect in range(nsystems):
    
    ax1.add_patch(patches.Rectangle( (amin,allypositions[iselect][0]-0.5*patchwidth), amax, patchwidth, color='gray', alpha=0.2))
    
    ax1.scatter(allsystems[iselect][:,finalcoldict['a']],allypositions[iselect][:],s=allsizes[iselect][:], color=allcolors[iselect][:])

    ax1.errorbar(allsystems[iselect][:,finalcoldict['a']], allypositions[iselect], xerr = allsystems[iselect][:,finalcoldict['a']]*allsystems[iselect][:,finalcoldict['e']], linestyle='None', color='black')
    

ax1.set_xlim(1.0e-1,amax)
ax1.set_xscale('log')
ax1.set_ylim(0,yposition[0]+patchwidth)
ax1.set_xlabel('a (AU)')
ax1.get_yaxis().set_ticks([])
plt.show()
fig1.savefig('portrait.png')

