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


nsystems = 100 # Number of systems to plot
ntrymax = 100000 # Number of samplings to try
systemspacing = 1000
patchwidth = 0.7*systemspacing
figheight = nsystems
figwidth = 0.8*nsystems

plt.rcParams['font.size'] = 18*np.int(nsystems)/10
plt.rcParams['patch.linewidth'] = 2
plt.rc('xtick', labelsize=14*np.int(nsystems)/10)
plt.rc('ytick', labelsize=14*np.int(nsystems)/10)

markerscale = 3*nsystems/10
minmarkersize = 50*nsystems/10
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

minsemimaj=[]

ntries = 0
iselect = 0
while isystem < nsystems:

    #iselect = np.random.randint(0,high=imax)
    iselect = iselect+1
    ntries = ntries+1
    if(ntries>ntrymax): break
    # Find first instance of iselect

    output = finaldata[finaldata[:,finalcoldict['istar']]==iselect,:]

    nplanets = output.shape[0]
    
    # Skip empty systems
    if nplanets<1: continue

    # Track minimum semimajor axis for sorting plots
    
    minsemimaj.append(np.amin(output[:,finalcoldict['a']]))


    
    # Check system hasn't been picked before
    if iselect in iselected: continue
    
    print 'Selected system ', iselect
    iselected.append(iselect)
    allsystems.append(output)
    
    isystem = isystem+1

    # Marker sizes
    sizes = markerscale*output[:,finalcoldict['mass']]

    sizes[sizes<minmarkersize] = minmarkersize

    allsizes.append(sizes)

    colors = []
   
   # print "Classifying Fragments"

    for i in range(nplanets):
        colors.append(get_fragment_colour(output,i))
    
    allcolors.append(colors)

# Now sort arrays in order of ascending minimum semimajor axis

neworder = np.flip(np.argsort(minsemimaj))
print neworder

print [minsemimaj[i] for i in neworder]

#print [allsystems[i][:,finalcoldict['a']] for i in range(nsystems)]
#allsystems.sort(key=lambda x: x[finalcoldict['a']])
#print [allsystems[i][:,finalcoldict['a']] for i in range(nsystems)]


#allsystems = [allsystems[i] for i in neworder]
#allypositions = [allypositions[i] for i in neworder]
#allsizes = [allsizes[i] for i in neworder]
#allcolors = [allcolors[i] for i in neworder]


for j in range(nsystems):
    
    
    i = neworder[j]
    nplanets = len(allsystems[i][:,finalcoldict['a']])
    #i=j
    
    # Generate y positions for each system
    ypositions = np.zeros(nplanets)
    ypositions[:] = j*systemspacing
    
    print i, nplanets, ypositions[0], allsystems[i][:,finalcoldict['a']], minsemimaj[i]
    
    ax1.add_patch(patches.Rectangle( (amin,ypositions[0]-0.5*patchwidth), amax, patchwidth, color='gray', alpha=0.2))
    
    ax1.scatter(allsystems[i][:,finalcoldict['a']],ypositions[:],s=allsizes[i][:], color=allcolors[i][:])

    ax1.errorbar(allsystems[i][:,finalcoldict['a']], ypositions[:], xerr = allsystems[i][:,finalcoldict['a']]*allsystems[i][:,finalcoldict['e']], linestyle='None', color='black')
    

ax1.set_xlim(1.0e-1,amax)
ax1.set_xscale('log')
#ax1.set_ylim(0,ypositions[0]+patchwidth)
ax1.set_xlabel('a (AU)')
ax1.get_yaxis().set_ticks([])
#plt.show()
fig1.savefig('portrait.png')

