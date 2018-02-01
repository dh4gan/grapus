# Written 09/03/2017 by dh4gan
# Compares the output from two runs of grapus v3.0 (nbody calculations)
# Makes scatter plots and histograms

import numpy as np
import matplotlib.pyplot as plt
import filefinder as ff
from os.path import splitext as ext
from io_grapus import ninitialcol,nfinalcol, initialcoldict,finalcoldict,finallabeldict, logcoldict


print ''
print '\t \t \t ------------------------'
print '\t \t \t Comparing Tidal Downsizing Runs v3.0 (N body data)'
print '\t \t \t ------------------------'
print ''

plt.rcParams['font.size'] = 18
plt.rcParams['patch.linewidth'] = 2
plt.rc('xtick', labelsize=14)
plt.rc('ytick', labelsize=14)

#prefix1 = 'control_nbodyon/total'
#prefix2 = 'feedback_nbodyon/total'

#run1name = 'No Core Feedback'
#run2name = 'Core Feedback'

#run1name = r'$C_{\rm space} = 3$'
#run2name = r'$C_{\rm space} = 10$'
#prefix1 = 'control_nbodyon/total'
#prefix2 = 'control_nbodyon_maxsep10/total'

#run1name = 'Initial $e=0$'
#run2name = 'Initial $e$ from Hall et al (2017)'
#prefix1='control_nbodyon/total'
#prefix2='control_initialecc/total'

#run1name = 'Control'
#run2name = r'$C_{\rm mig}=0.1$'
#prefix1='control_nbodyon/total'
#prefix2='Cmig_1e-1_nbodyon/total'

#run1name = 'Control'
#run2name = r'no Gap Formation'
#prefix1 = 'control_nbodyon/total'
#prefix2 = 'no_gaps/total'

run1name = 'Control'
run2name = r'$C_{\rm mig}=0.1$, no Gap Formation'
prefix1 = 'control_nbodyon/total'
prefix2 = 'Cmig_1e-1_nogaps/total'

# Read in initial file

planetstring = finallabeldict['mass']
astring = finallabeldict['a']
corestring = finallabeldict['mcore']

initialfile1 = prefix1+'.initial'
initialfile2 = prefix2+'.initial'

print 'Reading ',initialfile1

initialdata1 = np.genfromtxt(initialfile1)
ninitialrow1 = initialdata1.size/ninitialcol

print 'There are ',ninitialrow1, ' rows'

initialdata1.reshape(ninitialrow1,ninitialcol)


print 'Reading ',initialfile2

initialdata2 = np.genfromtxt(initialfile2)
ninitialrow2 = initialdata2.size/ninitialcol

print 'There are ',ninitialrow2, ' rows'

initialdata2.reshape(ninitialrow2,ninitialcol)



acol = initialcoldict['a']
mcol = initialcoldict['mass']
sedcol = initialcoldict['tsed']

initialdata1 = initialdata1[initialdata1[:,acol]>0.1]
initialdata2 = initialdata2[initialdata2[:,acol]>0.1]

print 'Initial Data for both runs read'

# Plot initial data

print 'Plotting mass versus semimajor axis for initial embryo distribution'

# Plot histograms

print 'Plotting histograms of initial distribution'

# Find maximum mass, and location

maxmass1 = np.amax(initialdata1[:,mcol])
maxval1 = np.argmax(initialdata1[:,mcol])

maxmass2 = np.amax(initialdata2[:,mcol])
maxval2 = np.argmax(initialdata2[:,mcol])

# Mass

mass = plt.figure()
ax = mass.add_subplot(111)
ax.set_xlabel(planetstring)

ax.hist(initialdata1[:,mcol],bins=100, label=run1name,histtype='step',color='blue', normed=True)
ax.hist(initialdata2[:,mcol],bins=100, label=run2name,histtype='step',color='red', linestyle='dashed',normed=True)
ax.legend()
mass.savefig('tworuns_mhist_initial.png', format='png')

# Semi major axis

sma = plt.figure()
ax = sma.add_subplot(111)
ax.set_xlabel(astring)
ax.set_xlim(0,1000.0)
ax.hist(initialdata1[:,acol],bins=100,range=[0,1000.0], label=run1name)
ax.hist(initialdata2[:,acol],bins=100,range=[0,1000.0], label=run2name)
ax.legend()
sma.savefig('tworuns_ahist_initial.png',format='png')

# Read in final files

acol = finalcoldict['a']
ecol = finalcoldict['e']
inccol = finalcoldict['i']
mcol = finalcoldict['mass']
mcorecol = finalcoldict['mcore']

finalfile1=prefix1+'.final'

print 'Reading ',finalfile1

finaldata1 = np.genfromtxt(finalfile1)
nfinalrow1 = finaldata1.size/nfinalcol
finaldata1.reshape(nfinalrow1,nfinalcol)

print 'There are initially ',nfinalrow1, 'rows'

# Delete rows where semi major axis equal to the boundary value rin

finaldata1 = finaldata1.compress(finaldata1[:,acol]>0.11, axis=0)
nfinalrow1 = finaldata1.shape[0]

print 'After deletion at boundary, there are ',nfinalrow1, ' rows'

# Separate data into ejected and non-ejected bodies

ejectadata1 = finaldata1[finaldata1[:,ecol]>=1.0]
nejected1 = ejectadata1.shape[0]

finaldata1 = finaldata1[finaldata1[:,ecol]<1.0]
nbound1 = finaldata1.shape[0]

print 'There are ',nejected1, ' ejected bodies'
print 'There are ',nbound1, ' bodies still bound to their parent stars'


finalfile2=prefix2+'.final'

print 'Reading ',finalfile2

finaldata2 = np.genfromtxt(finalfile2)
nfinalrow2 = finaldata2.size/nfinalcol
finaldata2.reshape(nfinalrow2,nfinalcol)

print 'There are initially ',nfinalrow2, 'rows'

# Delete rows where semi major axis equal to the boundary value rin

finaldata2 = finaldata2.compress(finaldata2[:,acol]>0.11, axis=0)
nfinalrow2 = finaldata2.shape[0]

print 'After deletion at boundary, there are ',nfinalrow2, ' rows'

# Separate data into ejected and non-ejected bodies

ejectadata2 = finaldata2[finaldata2[:,ecol]>=1.0]
nejected2 = ejectadata2.shape[0]

finaldata2 = finaldata2[finaldata2[:,ecol]<1.0]
nbound2 = finaldata2.shape[0]

print 'There are ',nejected2, ' ejected bodies'
print 'There are ',nbound2, ' bodies still bound to their parent stars'


print 'Data read'


if nejected1*nejected2 > 0.0:

    # Plot ejecta mass histogram

    ejmass = plt.figure()
    ax = ejmass.add_subplot(111)
    ax.set_xlabel(r'Ejecta Mass ($M_{\rm Jup}$)',fontsize=22)
    ax.set_ylabel('Relative Frequency',fontsize=22)
    ax.hist(ejectadata1[:,finalcoldict['mass']], bins=100, label=run1name, histtype = 'step', edgecolor='blue', normed=True)
    ax.hist(ejectadata2[:,finalcoldict['mass']], bins=100, label=run2name, histtype = 'step', edgecolor='green', normed=True,linestyle='dashed')
    ax.legend()
    ejmass.savefig('tworuns_ejectamass_hist.png',format= 'png')


    # Find ejecta hyperbolic velocity - needs mass of parent star

    print 'Finding stellar masses for run 1 : reading log file'

    logfile1 = prefix1+'.log'

    logdata1 = np.genfromtxt(logfile1)

    print 'Log file ', logfile1, ' read'

    ejectavelocity1 = np.zeros(nejected1)

    print 'Calculating ejecta velocities (at infinity)'

    for ieject in range(nejected1):
        # Find stellar mass for this body
        starindex = np.where(logdata1[:,logcoldict['istar']]==ejectadata1[ieject,finalcoldict['istar']])
        
        # Calculate ejecta velocity at infinity (converting to km s^-1)
        ejectavelocity1[ieject] = 29.84*np.sqrt(logdata1[starindex,logcoldict['mstar']]/finaldata1[ieject,finalcoldict['a']])
    
        #print logdata[starindex,logcoldict['mstar']], finaldata[ieject,finalcoldict['a']], ejectavelocity[ieject]

    print 'Finding stellar masses for run 2 : reading log file'

    logfile2 = prefix2+'.log'

    logdata2 = np.genfromtxt(logfile2)

    print 'Log file ', logfile2, ' read'

    ejectavelocity2 = np.zeros(nejected2)

    print 'Calculating ejecta velocities (at infinity)'

    for ieject in range(nejected2):
        # Find stellar mass for this body
        starindex = np.where(logdata2[:,logcoldict['istar']]==ejectadata2[ieject,finalcoldict['istar']])
   
        # Calculate ejecta velocity at infinity (converting to km s^-1)
        ejectavelocity2[ieject] = 29.84*np.sqrt(logdata2[starindex,logcoldict['mstar']]/finaldata2[ieject,finalcoldict['a']])
    
        #print logdata[starindex,logcoldict['mstar']], finaldata[ieject,finalcoldict['a']], ejectavelocity[ieject]



    print 'Plotting'
    # Plot distribution of ejecta velocity

    ejvel = plt.figure()
    ax = ejvel.add_subplot(111)
    ax.set_xlim(0,20)
    ax.set_xlabel('Ejecta Velocity ($km\, s^{-1}$)',fontsize=22)
    ax.set_ylabel('Relative Frequency',fontsize=22)
    ax.hist(ejectavelocity1, bins=100, range=[0,30], label=run1name, histtype = 'step', edgecolor='blue',normed=True)
    ax.hist(ejectavelocity2, bins=100, range=[0,30], label=run2name, histtype = 'step', edgecolor='red',normed=True,linestyle='dashed')
    ax.legend()
    ejvel.savefig('tworuns_ejectavel_hist.png',format= 'png')

# Plot final data

ncores1 = np.count_nonzero(finaldata1[:,mcorecol])
ncores2 = np.count_nonzero(finaldata2[:,mcorecol])
print 'There are ', ncores1, ' cores in run 1'
print 'There are ', ncores2, ' cores in run 2'
    
# Plot histograms

print 'Plotting histograms of final distribution'

# Mass

mass = plt.figure()
ax = mass.add_subplot(111)
ax.set_xlabel(finallabeldict['mass'],fontsize=22)
ax.set_ylabel('Relative Frequency',fontsize=22)
ax.hist(finaldata1[:,mcol],bins=100,label=run1name,histtype='step',color='blue', normed=True)
ax.hist(finaldata2[:,mcol],bins=100,label=run2name,histtype='step',color='red', linestyle='dashed',normed=True)
ax.legend()
mass.savefig('tworuns_mhist_final.png', format='png')

# Semi major axis

sma = plt.figure()
ax = sma.add_subplot(111)
ax.set_xlabel(finallabeldict['a'],fontsize=22)
ax.set_ylabel('Relative Frequency',fontsize=22)
ax.hist(finaldata1[:,acol],bins=100,range=[0,500],label=run1name, histtype='step',color='blue', normed=True)
ax.hist(finaldata2[:,acol],bins=100,range=[0,500],label=run2name, histtype='step',color='#cc9900', linestyle='dashed', normed=True)
#ax.hist(finaldata2[:,acol],bins=100,range=[0,500],label=run2name, histtype='step',color='red', linestyle='dashed', normed=True)
ax.legend()
sma.savefig('tworuns_ahist_final.png',format='png')

# Eccentricity

ecc = plt.figure()
ax = ecc.add_subplot(111)
ax.set_xlabel(finallabeldict['e'],fontsize=22)
ax.set_ylabel('Relative Frequency',fontsize=22)
ax.hist(finaldata1[:,ecol],bins=100, label=run1name, histtype='step',color='blue',normed=True)
ax.hist(finaldata2[:,ecol],bins=100, label=run2name, histtype='step',color='#cc00cc', linestyle='dashed',normed=True)
#ax.hist(finaldata2[:,ecol],bins=100, label=run2name, histtype='step',color='red', linestyle='dashed',normed=True)
ax.legend()
ecc.savefig('tworuns_ehist_final.png',format='png')

# Inclination

inc = plt.figure()
ax = inc.add_subplot(111)
ax.set_xlabel(finallabeldict['i'])
ax.set_ylabel('Relative Frequency')
ax.hist(finaldata1[:,inccol],bins=100, label=run1name, histtype='step',color='blue')
ax.hist(finaldata2[:,inccol],bins=100, label=run2name, histtype='step',color='red', linestyle='dashed')
inc.savefig('tworuns_inchist_final.png',format='png')

# Core mass

if(ncores1+ncores2>0):
    coredata1 = finaldata1[:,mcorecol][finaldata1[:,mcorecol]>0.0]
    coredata2 = finaldata2[:,mcorecol][finaldata2[:,mcorecol]>0.0]
    
    coreembryos1 = finaldata1[:,mcol][finaldata1[:,mcorecol]>0.0]
    
    

    core = plt.figure()
    ax = core.add_subplot(111)
    ax.set_xlabel(corestring, fontsize=22)
    ax.set_ylabel('N',fontsize=22)
    ax.hist(coredata1,bins=100,label=run1name, histtype='step',color='blue')
    ax.hist(coredata2,bins=100,label=run2name, histtype='step',color='#6f6867')
    ax.legend(loc='upper left')
    core.savefig('tworuns_mcorehist_final.png',format='png')

print 'Plotting comparison histograms'


# Bar chart of types

print "Classifying Fragments"

names = ('Icy Cores', 'Rocky Cores','No Cores', 'Brown Dwarfs', 'Planetesimals')

types1 = np.zeros(5)
types2 = np.zeros(5)

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


for i in range(len(finaldata1[:,1])):

    # Convert data into binary code
    binary = 0
    for j in range(1,7):
        binary += int(finaldata1[i,j])*2**(j-1)    
    
    # Check for icy core
    #print binary,np.sum(np.argwhere(icy==binary))>0,np.sum(np.argwhere(rocky==binary))>0,np.sum(np.argwhere(nocore==binary))>0,np.sum(np.argwhere(planetesimal==binary))>0 
    
    # Icy Core
    if(np.sum(np.argwhere(icy==binary)>0)): types1[0] +=1
    
    # Rocky Core
    if(np.sum(np.argwhere(rocky==binary)>0)): types1[1] +=1
    
    # No Core /BD
    if(np.sum(np.argwhere(nocore==binary)>0)):
        # BD 
        if finaldata1[i,finalcoldict['mass']] >13.0 and finaldata1[i,3]==1:
            types1[3]+=1
        # No core
        else:
            types1[2] +=1        
    
    # Planetesimal Belt    
    if(np.sum(np.argwhere(planetesimal==binary)>0)): types1[4] +=1

    # Keep track of unclassified objects    
    if(not(np.sum(np.argwhere(alltypes==binary)>0))):
        if finaldata1[i,finalcoldict['mass']] > 13.0: 
            types1[3]+=1
        else:
            counter +=1        
            print 'Unclassified', binary, finaldata1[i,1:7], finaldata1[i,finalcoldict['mass']],finaldata1[i,finalcoldict['a']]
    
            
types1[:] = 100.0*types1[:]/len(finaldata1[:,1])

for i in range(len(types1)):
    print names[i],types1[i]
    
# Plot histogram

print 'There are ',counter, ' unclassified objects'

print 'Plotting bar chart'



fate=plt.figure()
ax=fate.add_subplot(111)
ax.set_xlabel('Relative Frequency (%)')
ax.set_ylabel('Type')
ax.set_ylim(0,6)

x = [0.5,1.5,2.5,3.5,4.5,5.5]



print 'The destruction rate is ', (float(nfinalrow1)/float(ninitialrow1))*100.0 , ' %'
print 'The ejection rate is ',(float(nejected1)/float(nfinalrow1))*100.0, ' %'

 
# Set up string tick labels
tick_locs = (x)
tick_lbls = names
plt.yticks(tick_locs, tick_lbls, rotation = 45)


ax.barh(x[0],types1[0],color='blue')
ax.barh(x[1],types1[1],color = 'red')
ax.barh(x[2],types1[2],color='green')
ax.barh(x[3],types1[3],color='brown')
ax.barh(x[4],types1[4],color='black')
fate.savefig("types.png", format='png')


