#
#
# Python module for grapus output data
# Helpful global variables and functions 
# for reading, writing data
#
#

import numpy as np

# Some helpful units
mearth_in_mjup = 1.0/317.8 # 1 Earth mass in Jupiter masses



# Initial File Data

ninitialcol = 8
initialcols = range(ninitialcol)
initialkeys = ['a', 'mass', 'radius', 'T0', 'scrit', 'tcool','tgrow','tsed']
initiallabels = ['Semimajor Axis (AU)', 'Initial Fragment Mass ($M_{Jup}$)', r'Initial Radius ($R_{Jup}$)', 'T0','scrit', 'tcool','tgrow','tsed']

initialcoldict =dict(zip(initialkeys,initialcols))
initiallabeldict =dict(zip(initialkeys,initiallabels))
initialfmt = '%.8e %.8e %.8e %.8e %.8e %.8e %.8e %.8e'


# Final (and snapshot) Data

nfinalcol=16

finalcols = range(nfinalcol)
finalkeys = ['istar', 'iembryo','imelt', 'ivap','idiss','igrown','iself', 'ijeans', 'itidal', 'a', 'e', 'i', 'mass', 'radius', 'mcore', 'rcore']
finallabels = ['Star Index', 'Embryo Index','imelt', 'ivap','idiss','igrown','iself', 'ijeans', 'itidal', 'Semimajor Axis (AU)', 'Eccentricity', 'Inclination', 'Final Fragment Mass ($M_{Jup}$) ', 'Final Fragment Radius ($R_{Jup}$)', 'Core Mass ($M_{\oplus}$)', 'Core Radius ($R_{\oplus}$)']

finalcoldict = dict(zip(finalkeys,finalcols))
finallabeldict = dict(zip(finalkeys,finallabels))

finalfmt = '%i %i %i %i %i %i %i %i %i %.8e %.8e %.8e %.8e %.8e %.8e %.8e'

nfinalcol = 15 # FOR USE ON OLD DATA
finalcoldict['a']= finalcoldict['a']-1
finalcoldict['e']= finalcoldict['e']-1
finalcoldict['mass']= finalcoldict['mass']-1
finalcoldict['mcore']= finalcoldict['mcore']-1

# Log Data

nlogcol=7
logcols = range(nlogcol)
logkeys = ['istar', 'mstar', 'mdisc', 'q_disc', 'rout', 'rfrag', 'nembryo']
loglabels = ['istar', 'mstar', 'mdisc', 'q_disc', 'rout', 'rfrag', 'nembryo']

logcoldict = dict(zip(logkeys,logcols))
loglabeldict = dict(zip(logkeys,loglabels))

logfmt = '%i %.8e % .8e %.8e %.8e %.8e %i'


#
# Classifying the final state of disc fragments given their evolutionary state
#


# The below numbers are either zero or one
# (imelt ivap idiss igrown iself ijeans itidal)
# Each possible outcome represented by a set of binary numbers
# (listed below)
#

typelabels = ('Gas Giant, \nIcy Core', 'Gas Giant, \nRocky Core','Gas Giant, \nNo Core', 'Terrestrial Planet', 'Brown Dwarfs', 'Planetesimals')

#
# Objects with icy cores
# (0 0 0 1 1 1 *)  in binary = (56,120)
#

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


# Colour scheme for plotting system portraits

BDcolour = '#663300'  # Brown Dwarfs
giantcolour = 'red'        # Giant Planets
terrestrialcolour = '#0099ff'   # Terrestrial Planets

colourdict = {'brown dwarf':BDcolour,'giant':giantcolour, 'terrestrial':terrestrialcolour, '':'grey'}


def read_finaldata(finalfile):


    print 'Reading final population data from ',finalfile

    finaldata = np.genfromtxt(finalfile)
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

    finaldata = finaldata.compress(finaldata[:,finalcoldict['a']]>0.11, axis=0)

    # Delete low mass bodies without cores

    blob = (finaldata[:,finalcoldict['mass']]<1.0e-2) & (finaldata[:,finalcoldict['mcore']]==0.0)

    finaldata = finaldata[np.logical_not(blob)]

    nfinalrow = finaldata.shape[0]

    print 'After deletion at inner system boundary, there are ',nfinalrow, ' rows'


    # Separate data into ejected and non-ejected bodies

    ejectadata = finaldata[finaldata[:,ecol]>=1.0]
    nejected = ejectadata.shape[0]

    finaldata = finaldata[finaldata[:,ecol]<1.0]
    nbound = finaldata.shape[0]

    print 'There are ',nejected, ' ejected bodies'
    print 'There are ',nbound, ' bodies still bound to their parent stars'

    return finaldata, ejectadata, nbound,nejected


def return_system(finaldata,istar):
    '''Returns collected data for planetary system istar'''

    systemdata = finaldata[finaldata[:,finalcoldict['istar']]==istar,:]

    masses = systemdata[:,finalcoldict['mass']]
    semimaj = systemdata[:,finalcoldict['a']]
    ecc = systemdata[:,finalcoldict['e']]

    colours = []

    for i in range(len(systemdata[:,finalcoldict['istar']])):
        colours.append(get_fragment_colour(systemdata,i))


    return masses, semimaj, ecc, colours


def aggregate_systems(finaldata, systemlist=None):
    '''Collates mass, semimajor axis and eccentricity data by system

    finaldata - output data from grapus
    systemlist - a list of istar values, identifying systems to aggregate'''

    nsystems = int(np.amax(finaldata[:,finalcoldict['istar']]))

    print 'Aggregating ', nsystems, ' systems'

    allmasses = []
    allsemimaj = []
    allecc = []
    allcolours = []

    if(systemlist==None):
        syslist = range(nsystems)
    else:
        syslist = systemlist

    first = True
    for i in syslist:

        masses,semimaj,ecc,colours = return_system(finaldata,i)

        allmasses.append(masses)
        allsemimaj.append(semimaj)
        allecc.append(ecc)
        allcolours.append(colours)

    allmasses.pop(0)
    allsemimaj.pop(0)
    allecc.pop(0)

    return allmasses,allsemimaj,allecc



def classify_fragment(fragdata,i):
    '''Classifies a fragment into brown dwarf, gas giant or terrestrial'''
    # Convert data into binary code

    fragmentclass = ''

    binary = 0
    for j in range(1,7):
        binary += int(fragdata[i,j])*2**(j-1)    
    
    
    #
    # Icy Core
    #
    if(np.sum(np.argwhere(icy==binary)>0)): fragmentclass='giant, icy core'
    
    #
    # Rocky Core
    #
    if(np.sum(np.argwhere(rocky==binary)>0)):
        
        # Is this gas giant or terrestrial planet? determine by mcore/membryo
        mratio = fragdata[i,finalcoldict['mcore']]*mearth_in_mjup/fragdata[i,finalcoldict['mass']]
         
         #print mratio, fragdata[i,finalcoldict['mcore']], fragdata[i,finalcoldict['mass']]
        if(mratio < 0.5):
            fragmentclass = 'giant, rocky core'
        else:
            fragmentclass = 'terrestrial'
    

    #
    # No Core /BD
    #
    if(np.sum(np.argwhere(nocore==binary)>0)):
            # BD 
            if (fragdata[i,finalcoldict['mass']] >13.0 and fragdata[i,finalcoldict['idiss']]==1):
                fragmentclass = 'brown dwarf'
                # No core
            else:
                fragmentclass = 'giant, no core'

    if(fragmentclass=='terrestrial'):
        print fragdata[i,1:7], binary, fragmentclass
    return fragmentclass


def get_fragment_colour(fragdata,i):
    
    fragmentclass = classify_fragment(fragdata,i)

    if('giant' in fragmentclass):
        fragmentclass = 'giant'

    return colourdict[fragmentclass]
    
