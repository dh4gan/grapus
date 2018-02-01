# Written 09/01/2017 by dh4gan
# Plots the output from grapus v3.0 (nbody calculations)


import numpy as np
import matplotlib.pyplot as plt
import filefinder as ff
from os.path import splitext as ext
from collections import Counter
from io_grapus import ninitialcol,nfinalcol, initialcoldict,finalcoldict,finallabeldict, logcoldict
from colormaps import viridis


# Python script reads in grapus data and does scatter plots and histograms


plt.rcParams['font.size'] = 18
plt.rcParams['patch.linewidth'] = 2
plt.rc('xtick', labelsize=14)
plt.rc('ytick', labelsize=14)


print ''
print '\t \t \t ------------------------'
print '\t \t \t Tidal Downsizing Plotter v3.0 (N body data)'
print '\t \t \t ------------------------'
print ''

# Use filefinder to find file prefixes of interest to the user

prefix = ff.find_sorted_local_input_files('*.initial')

# Strip the file extension
prefix = ext(prefix)[0]

# Read in initial file

planetstring = finallabeldict['mass']
astring = finallabeldict['a']
corestring = finallabeldict['mcore']

initialfile=prefix+'.initial'

print 'Reading ',initialfile

initialdata = np.genfromtxt(initialfile)
ninitialrow = initialdata.size/ninitialcol

print 'There are ',ninitialrow, ' rows'

initialdata.reshape(ninitialrow,ninitialcol)


acol = initialcoldict['a']
mcol = initialcoldict['mass']
sedcol = initialcoldict['tsed']

initialdata = initialdata[initialdata[:,acol]>0.1]

print 'Data read'

# Plot initial data

print 'Plotting mass versus semimajor axis for initial embryo distribution'

# Mass vs Semi Major Axis

mvsa=plt.figure()
ax = mvsa.add_subplot(111)
ax.set_xlim(1.0e-1,1.0e2)
#ax.set_ylim(1.0e0,1e3)
ax.set_xlabel(astring)
ax.set_ylabel(planetstring)
ax.set_xscale('log')
ax.set_yscale('log')

scat = ax.scatter(initialdata[:,acol],initialdata[:,mcol],c=initialdata[:,sedcol])
colors= mvsa.colorbar(scat,cax=None)
colors.set_label("Sedimentation Time (yr)")

mvsa.savefig('m_vs_a_initial.png', format='png')

print 'Saved to m_vs_a_initial.png'

# Plot histograms

print 'Plotting histograms of initial distribution'

# Find maximum mass, and location

maxmass = np.amax(initialdata[:,mcol])
maxval = np.argmax(initialdata[:,mcol])

# Mass

mass = plt.figure()
ax = mass.add_subplot(111)
ax.set_xlabel(planetstring)

histo = ax.hist(initialdata[:,mcol],bins=100)
mass.savefig('mhist_initial.png', format='png')

# Semi major axis

sma = plt.figure()
ax = sma.add_subplot(111)
ax.set_xlabel(astring)
histo = ax.hist(initialdata[:,acol],bins=100)
sma.savefig('ahist_initial.png',format='png')

# Read in final file
finalfile=prefix+'.final'

print 'Reading ',finalfile

finaldata = np.genfromtxt(finalfile)
nfinalrow = finaldata.size/nfinalcol
finaldata.reshape(nfinalrow,nfinalcol)

print 'There are initially ',nfinalrow, 'rows'

# Delete rows where semi major axis equal to the boundary value rin

acol = finalcoldict['a']
ecol = finalcoldict['e']
inccol = finalcoldict['i']
mcol = finalcoldict['mass']
mcorecol = finalcoldict['mcore']

finaldata = finaldata.compress(finaldata[:,acol]>0.11, axis=0)
# Delete low mass bodies without cores

blob = (finaldata[:,mcol]<1.0e-2) & (finaldata[:,mcorecol]==0.0)

finaldata = finaldata[np.logical_not(blob)]
#finaldata = finaldata.compress(not(blob), axis=0)

nfinalrow = finaldata.shape[0]

print 'After deletion at boundary, there are ',nfinalrow, ' rows'

print 'Data read'

# Separate data into ejected and non-ejected bodies

ejectadata = finaldata[finaldata[:,ecol]>=1.0]
nejected = ejectadata.shape[0]

finaldata = finaldata[finaldata[:,ecol]<1.0]
nbound = finaldata.shape[0]

print 'There are ',nejected, ' ejected bodies'
print 'There are ',nbound, ' bodies still bound to their parent stars'

logfile = prefix+'.log'

logdata = np.genfromtxt(logfile)

print 'Log file ', logfile, ' read'

if (nejected>0):
    # Plot ejecta mass histogram

    ejmass = plt.figure()
    ax = ejmass.add_subplot(111)
    ax.set_xlabel(r'Mass ($M_{\rm Jup}$)')
    ax.set_ylabel('N')
    ax.hist(ejectadata[:,finalcoldict['mass']], bins=100)
    ejmass.savefig('ejectamass_hist.png',format= 'png')


    # Find ejecta hyperbolic velocity - needs mass of parent star

    print 'Finding stellar masses: reading log file'

    

    ejectavelocity = np.zeros(nejected)

    print 'Calculating ejecta velocities (at infinity)'

    for ieject in range(nejected):
        # Find stellar mass for this body
        starindex = np.where(logdata[:,logcoldict['istar']]==ejectadata[ieject,finalcoldict['istar']])
   
        # Calculate ejecta velocity at infinity (converting to km s^-1)
        ejectavelocity[ieject] = 29.84*np.sqrt(logdata[starindex,logcoldict['mstar']]/finaldata[ieject,finalcoldict['a']])
    
        #print logdata[starindex,logcoldict['mstar']], finaldata[ieject,finalcoldict['a']], ejectavelocity[ieject]

    print 'Plotting'
    # Plot distribution of ejecta velocity

    ejvel = plt.figure()
    ax = ejvel.add_subplot(111)
    ax.set_xlabel('Ejecta Velocity ($km\, s^{-1}$)')
    ax.set_ylabel('N')
    ax.hist(ejectavelocity, bins=100)
    ejvel.savefig('ejectavel_hist.png',format= 'png')

# Plot final data

# Count system multiplicity

print 'Obtaining multiplicity data'

multiplicitycounts = Counter(finaldata[:,finalcoldict['istar']])
initialmultiplicity = logdata[:,logcoldict['nembryo']]

maxstar = int(logdata[-1,logcoldict['istar']])

val = []
for i in range(maxstar):    
    val.append(multiplicitycounts[float(i)])    
    
counts = plt.figure()
ax = counts.add_subplot(111)
ax.set_xlabel('Multiplicity', fontsize=20)
ax.set_ylabel('Relative Frequency',fontsize=20)
ax.set_ylim(0,7)
ax.hist(val, normed=True, bins=20,label='Final',histtype='step',edgecolor='r')
ax.hist(initialmultiplicity, normed=True,bins=20,histtype='step',edgecolor = 'b',linestyle='dotted',label='Initial')
ax.legend()
counts.savefig('multiplicity.png')

# Begin with m_vs_a

print 'Plotting mass versus semimajor axis for final (bound) distribution'

mvsa=plt.figure()

ax=mvsa.add_subplot(111)
ax.set_xlim(1.0e-1,1.0e5)
ax.set_ylim(1.0e-6,1e3)
#ax.set_xlim(1.0e-1,1.0e3)
#ax.set_ylim(1.0e-1,1e3)
ax.set_xlabel(finallabeldict['a'])
ax.set_ylabel(finallabeldict['mass'])
ax.set_xscale('log')
ax.set_yscale('log')

# Take size data and scale correctly

colourticks=np.linspace(0.0,10.0<1.0,endpoint=True)

scat = ax.scatter(finaldata[:,acol],finaldata[:,mcol],s=(finaldata[:,mcorecol]*10)+10, c=finaldata[:,mcorecol],cmap=viridis)

colours= mvsa.colorbar(scat,cax=None)
colours.set_label(corestring)

mvsa.savefig('m_vs_a_final.png', format='png')

# Now plot e vs a

print 'Plotting eccentricity versus semimajor axis for final (bound) distribution'

mvsa=plt.figure()

ax=mvsa.add_subplot(111)
#ax.set_xlim(1.0e-3,1.0e3)
#ax.set_ylim(1.0e-4,1.0e0)
ax.set_xlabel(astring)
ax.set_ylabel('Eccentricity')
ax.set_xscale('log')
ax.set_yscale('log')

# Take size data and scale correctly

colourticks=np.linspace(0.0,10.0<1.0,endpoint=True)

#scat = ax.scatter(finaldata[:,7],finaldata[:,8], s=finaldata[:,9],c=finaldata[:,11], alpha=0.5)
scat = ax.scatter(finaldata[:,acol],finaldata[:,ecol],s=(finaldata[:,mcorecol]*10)+10, c=finaldata[:,mcorecol])
#scat2 = ax.scatter(exodata[:,0], exodata[:,1],marker='o', edgecolor = 'black', facecolor = 'none')
colours= mvsa.colorbar(scat,cax=None)
colours.set_label(corestring)

mvsa.savefig('e_vs_a_final.png', format='png')

ncores = np.count_nonzero(finaldata[:,mcorecol])
print 'There are ', ncores, ' cores'

# m vs e

mvse=plt.figure()

ax=mvse.add_subplot(111)
#ax.set_xlim(1.0e-3,1.0e3)
#ax.set_ylim(1.0e-4,1.0e0)
ax.set_xlabel(finallabeldict['mass'])
ax.set_ylabel('Eccentricity')
ax.set_xscale('log')
ax.set_yscale('log')

# Take size data and scale correctly

colourticks=np.linspace(0.0,10.0<1.0,endpoint=True)

#scat = ax.scatter(finaldata[:,7],finaldata[:,8], s=finaldata[:,9],c=finaldata[:,11], alpha=0.5)
scat = ax.scatter(finaldata[:,mcol],finaldata[:,ecol],s=(finaldata[:,mcorecol]*10)+10, c=finaldata[:,mcorecol])
#scat2 = ax.scatter(exodata[:,0], exodata[:,1],marker='o', edgecolor = 'black', facecolor = 'none')
colours= mvsa.colorbar(scat,cax=None)
colours.set_label(corestring)

mvse.savefig('m_vs_e_final.png', format='png')


# Plot histograms

print 'Plotting histograms of final distribution'

# Mass

mass = plt.figure()
ax = mass.add_subplot(111)
ax.set_xlabel(finallabeldict['mass'])
ax.set_ylabel('N')
histo = ax.hist(finaldata[:,mcol],bins=100)
mass.savefig('mhist_final.png', format='png')

# Semi major axis

sma = plt.figure()
ax = sma.add_subplot(111)
ax.set_xlabel(finallabeldict['a'])
ax.set_ylabel('N')
histo = ax.hist(finaldata[:,acol],bins=100)
sma.savefig('ahist_final.png',format='png')

# Eccentricity

ecc = plt.figure()
ax = ecc.add_subplot(111)
ax.set_xlabel(finallabeldict['e'])
ax.set_ylabel('N')
histo = ax.hist(finaldata[:,ecol][finaldata[:,ecol]<1.0],bins=100)
ecc.savefig('ehist_final.png',format='png')

# Inclination

inc = plt.figure()
ax = inc.add_subplot(111)
ax.set_xlabel(finallabeldict['i'])
ax.set_ylabel('N')
ax.hist(finaldata[:,inccol],bins=100)
inc.savefig('inchist_final.png',format='png')


# Core mass

if(ncores>0):
    
    acore = finaldata[:,acol][finaldata[:,mcorecol]>0.0]
    coremasses = finaldata[:,mcorecol][finaldata[:,mcorecol]>0.0]
    coreembryos = finaldata[:,mcol][finaldata[:,mcorecol]>0.0]
    
    corefraction = coremasses*0.003146/(coremasses+coreembryos)

    core = plt.figure()
    ax = core.add_subplot(111)
    ax.set_xlabel(corestring)
    ax.set_ylabel('N')
    histo = ax.hist(coremasses,bins=50)
    core.savefig('mcorehist_final.png',format='png')

    core = plt.figure()
    ax = core.add_subplot(111)
    ax.set_xlabel(corestring)
    ax.set_ylabel('Core to Object Mass Ratio')
    ax.set_yscale('log')
    ax.set_xscale('log')
    ax.set_xlim(1.0e-1,1.0e2)
    ax.set_ylim(1e-4,1e0)
    plot = ax.scatter(coremasses,coreembryos)
    core.savefig('m_vs_mcore_final.png', format='png')
    
    print 'Plotting mcore versus semimajor axis for final distribution'

    mvsa=plt.figure()

    ax=mvsa.add_subplot(111)
    #ax.set_xlim(1.0e-2,1.0e2)
    #ax.set_ylim(1.0e-1,1e2)
    ax.set_xlabel(astring)
    ax.set_ylabel(corestring)
    ax.set_xscale('log')
    ax.set_yscale('log')
    
    print np.amax(corefraction), np.amin(corefraction)

    scat = ax.scatter(acore,coremasses, c=corefraction, edgecolor='black', cmap='Blues')
    colours= mvsa.colorbar(scat,cax=None)
    colours.set_label('Core Mass Fraction')   
    mvsa.savefig('mcore_vs_a_final.png', format='png')


print 'Plotting comparison histograms'

mass = plt.figure()
ax = mass.add_subplot(111)
ax.set_xlabel(planetstring)
ax.set_ylabel('Relative Frequency')
#ax.hist([initialdata[:,initialcoldict['mass']], finaldata[:,mcol], ejectadata[:,mcol]],bins=100,stacked=True,normed=True,histtype='step', label=['Initial','Bound','Ejected'], color = ['b','r','g'])
histo = ax.hist(initialdata[:,initialcoldict['mass']],bins=100,histtype='step',edgecolor = 'b',linestyle='dotted',label='Initial',normed=False)
histo = ax.hist(finaldata[:,finalcoldict['mass']],bins=100,histtype='step',edgecolor='r',label='Bound',normed=False)
if(nejected>0): histo = ax.hist(ejectadata[:,finalcoldict['mass']],bins=100,histtype='step',edgecolor='g',linestyle='dashed',label='Ejected',normed=False)

ax.legend(loc='upper right')
mass.savefig('mhist_compare.png', format='png')

sma = plt.figure()
ax = sma.add_subplot(111)
ax.set_xlabel(astring)
ax.set_ylabel('Relative Frequency')
histo = ax.hist(initialdata[:,initialcoldict['a']],bins=100,histtype='step',edgecolor = 'b',label='Initial',normed=True)
histo = ax.hist(finaldata[:,finalcoldict['a']][finaldata[:,finalcoldict['a']]<1.0e3],bins=100,histtype='bar',color='r',label='Final',normed=True)
ax.legend(loc='upper right')
sma.savefig('ahist_compare.png', format='png')

print 'Plotting BD only histogram'

sma = plt.figure()
ax = sma.add_subplot(111)
ax.set_xlabel(astring)
ax.set_ylabel('Relative Frequency')
ax.set_xlim(0,100)
#histo = ax.hist(initialdata[:,0][initialdata[:,1]>13.0],bins=100,histtype='step',edgecolor = 'b',label='Initial')
histo = ax.hist(finaldata[:,finalcoldict['mass']][finaldata[:,finalcoldict['mass']]>13.0],bins=100,histtype='bar',color='r',normed=True)
#ax.legend(loc='upper right')
sma.savefig('ahist_BD.png', format='png')

# Bar chart of types

print "Classifying Fragments"

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


for i in range(len(finaldata[:,1])):

    # Convert data into binary code
    binary = 0
    for j in range(1,7):
        binary += int(finaldata[i,j])*2**(j-1)    
    
    # Check for icy core
    #print binary,np.sum(np.argwhere(icy==binary))>0,np.sum(np.argwhere(rocky==binary))>0,np.sum(np.argwhere(nocore==binary))>0,np.sum(np.argwhere(planetesimal==binary))>0 
    
    # Icy Core
    if(np.sum(np.argwhere(icy==binary)>0)): types[0] +=1
    
    # Rocky Core
    if(np.sum(np.argwhere(rocky==binary)>0)):
        
        # Is this gas giant or terrestrial planet? determine by mcore/membryo
        mratio = finaldata[i,mcorecol]*0.003146/finaldata[i,mcol]
         
        if(mratio < 0.5):
            types[1] +=1
        else:
            types[3]+=1
    
    # No Core /BD
    if(np.sum(np.argwhere(nocore==binary)>0)):
        # BD 
        if finaldata[i,finalcoldict['mass']] >13.0 and finaldata[i,3]==1:
            types[4]+=1
        # No core
        else:
            types[2] +=1        
    
    # Planetesimal Belt    
    if(np.sum(np.argwhere(planetesimal==binary)>0)): types[5] +=1

    # Keep track of unclassified objects    
    if(not(np.sum(np.argwhere(alltypes==binary)>0))):
        if finaldata[i,finalcoldict['mass']] > 13.0: 
            types[4]+=1
        else:
            counter +=1        
            print 'Unclassified', binary, finaldata[i,1:7], finaldata[i,finalcoldict['mass']],finaldata[i,finalcoldict['a']]
    
            
# Plot histogram

print 'There are ',counter, ' unclassified objects'

print 'Plotting bar chart'



fate=plt.figure()
ax=fate.add_subplot(111)
ax.set_xlabel('Relative Frequency (%)')
ax.set_ylabel('Type')
ax.set_ylim(0,6)

x = [0.5,1.5,2.5,3.5,4.5,5.5,6.5]

types[:] = 100.0*types[:]/len(finaldata[:,1])
print counter
for i in range(len(types)):
    print names[i],types[i]

print 'The destruction rate is ', (float(nfinalrow)/float(ninitialrow))*100.0 , ' %'
print 'The ejection rate is ',(float(nejected)/float(nfinalrow))*100.0, ' %'

 
# Set up string tick labels
tick_locs = (x)
tick_lbls = names
plt.yticks(tick_locs, tick_lbls, rotation = 60)


ax.barh(x[0],types[0],color='blue')
ax.barh(x[1],types[1],color = 'yellow')
ax.barh(x[2],types[2],color='red')
ax.barh(x[3],types[3],color='green')
ax.barh(x[4],types[4],color='brown')
ax.barh(x[5],types[5],color='black')
fate.savefig("types.png", format='png')

# Final step - write some useful numbers to a data file

print 'Writing global stats to globalstats.dat'

f = open ('globalstats.dat', 'w')

f.write('Total no. of fragments '+str(ninitialrow) + '\n')
f.write('Total no. of surviving fragments '+str(nfinalrow) + '\n')
f.write('Total no. of ejected fragments '+str(nejected) + '\n')
if(ncores>0):
    f.write('Total no. of cores '+str(len(coremasses)) + '\n')
else:
    f.write('Total no. of cores : 0  \n')

f.close()


