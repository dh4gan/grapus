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

finaldata, ejectadata, nbound,nejected = read_finaldata(finalfile)

isystem = 0

imax = int(np.amax(finaldata[:,finalcoldict['istar']]))

print 'Have ', imax, ' systems to choose from for portraits'

fig1 = plt.figure()
ax1 = fig1.add_subplot(111)

counter = 0

selected = []

ntries = 0
while isystem < nsystems:

    iselect = np.random.randint(0,high=imax)
    
    ntries = ntries+1
    if(ntries>ntrymax): break
    # Find first instance of iselect

    output = finaldata[finaldata[:,finalcoldict['istar']]==iselect,:]

    nplanets = output.shape[0]
    
    # Only plot systems with more than 3 bodies at a<10 AU
    if(nplanets<2 or np.amin(output[:,finalcoldict['a']])>10.0): continue
    
    
    # Check system hasn't been picked before
    if iselect in selected: continue
    
    print 'Selected system ', iselect
    selected.append(iselect)
    isystem = isystem+1

    yposition = np.zeros(nplanets)
    yposition[:]=isystem*systemspacing

    # Marker sizes
    sizes = markerscale*output[:,finalcoldict['mass']]

    sizes[sizes<minmarkersize] = minmarkersize

    colors = []
   

    print "Classifying Fragments"

    for i in range(nplanets):

        colors.append(get_fragment_colour(output,i))
                
        print output[:,finalcoldict['mass']], output[:,finalcoldict['mcore']],colors

    ax1.add_patch(patches.Rectangle( (amin,yposition[0]-0.5*patchwidth), amax, patchwidth, color='gray', alpha=0.2))
    ax1.scatter(output[:,finalcoldict['a']],yposition,s=sizes, color=colors)
    ax1.errorbar(output[:,finalcoldict['a']],yposition,xerr=output[:,finalcoldict['a']]*output[:,finalcoldict['e']], linestyle='None', color='black')
    

ax1.set_xlim(1.0e-1,amax)
ax1.set_xscale('log')
ax1.set_ylim(0,yposition[0]+patchwidth)
ax1.set_xlabel('a (AU)',fontsize=22)
ax1.get_yaxis().set_ticks([])
plt.show()
fig1.savefig('portrait.png')

