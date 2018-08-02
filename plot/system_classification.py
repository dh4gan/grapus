#
# Written 28/3/17 by dh4gan
# Uses unsupervised machine learning on systems
# generated from grapus v3.0
# to identify types
#

import matplotlib.pyplot as plt
import matplotlib.patches as patches
import filefinder as ff
import numpy as np
from io_grapus import read_finaldata,finalcoldict, aggregate_systems
from sklearn.cluster import KMeans
import sys


def add_variable_to_feature_vector(vector,bodies,nplanetmax,key):

    nbodies = bodies.shape[0]

    for i in range(nbodies):
        vector = np.append(vector,bodies[i,finalcoldict[key]])
    for i in range(nbodies,nplanetmax):
        vector = np.append(vector,0.0)

    return vector


def create_input_feature_vector(bodies,nplanetmax):
    '''generates a vector of data representing a planetary system for classification'''

    nbodies = bodies.shape[0]

    if(nbodies==0):
        print 'ZERO BODIES ', bodies
        
    vector = np.array([nbodies])

    max_semimaj = np.amax(bodies[:,finalcoldict['a']])
    max_mass = np.amax(bodies[:,finalcoldict['mass']])

    # Add semimajor axis, mass and eccentricity to feature vector

    vector = add_variable_to_feature_vector(vector,bodies,nplanetmax,'a')
    vector = add_variable_to_feature_vector(vector,bodies,nplanetmax,'mass')
    vector = add_variable_to_feature_vector(vector,bodies,nplanetmax,'e')
    
    return vector, max_semimaj, max_mass

print ''
print '\t \t \t ------------------------'
print '\t \t \t Tidal Downsizing System Classification (N body data)'
print '\t \t \t ------------------------'
print ''




plt.rcParams['font.size'] = 18
plt.rcParams['patch.linewidth'] = 2
plt.rc('xtick', labelsize=14)
plt.rc('ytick', labelsize=14)

nsystems = 10 # Number of systems to plot
nplanetmax = 3
ntrymax = 100000 # Number of samplings to try
systemspacing = 1000
patchwidth = 1000

markerscale = 3
minmarkersize = 50
amin = 0.1
amax=1.0e5

istarcol = finalcoldict['istar']
mcol = finalcoldict['mass']
acol = finalcoldict['a']

# Use filefinder to find files

finalfile = ff.find_sorted_local_input_files('*.final')

finaldata, ejectadata, nbound,nejected = read_finaldata(finalfile)

isystem = 0

imax = int(np.amax(finaldata[:,istarcol]))


print 'Have maximum of ', imax, ' systems to classify'

if nsystems>imax:
    nsystems=imax

fig1 = plt.figure()
ax1 = fig1.add_subplot(111)

print 'Building feature vectors for ',nsystems, ' systems'
# Now create vectors for each planetary system

systems = []

systemlist = range(10)
allmasses,allsemimaj,allecc = aggregate_systems(finaldata,systemlist)

print systemlist, allmasses


for isystem in range(1,nsystems):
    
    # Select planetary bodies
    bodies = finaldata[finaldata[:,istarcol]==isystem,:]

    #print output[:,mcorecol]*0.003146/output[:,mcol]
    nplanets = bodies.shape[0]    

    # Skip empty systems
    if(nplanets==0): continue

    # Skip systems which have more than nplanetmax planets
    if(nplanets>nplanetmax): continue
    
    systems.append(bodies)

    print 'SYSTEMS'
    print systems
    # Generate vector for this system
    output, max_a,max_m = create_input_feature_vector(bodies,nplanetmax)

    #print isystem,output

    # If this is the first system, generate the input matrix X
    if(isystem==1):
        X = output
        max_semimaj = np.array(max_a)
        max_mass = np.array(max_m)
        # Otherwise append system to pre-existing matrix
    else:
        X = np.vstack((X,output))
        max_semimaj = np.hstack((max_semimaj,max_a))
        max_mass = np.hstack((max_mass,max_m))



# Order systems list by maximum semimajor axis

print "ORGANISING"
print max_semimaj
print systems
systems = [x for _, x in sorted(zip(max_semimaj,systems))]

print "SORTED"
print systems


# Check input vector is OK

print "Feature vector built"

#print "X: ", X

nvectors = X.shape[0]
print nvectors

print  max_mass

# Now try k-means with various k-values

print nvectors, ' vectors constructed: classifying'

kvalues = [2,3,4,5]



colors = ['red','yellow','blue','green','purple', 'black']

for k in kvalues:
    print 'Testing k=',k

    fig1 = plt.figure(1)
    ax1 = fig1.add_subplot(111)
    kmeans = KMeans(n_clusters=k, n_init=10,random_state=0).fit(X)

    #print kmeans.cluster_centers_

    colorlist = []
    for i in range(nvectors):    
        colorlist.append(colors[kmeans.labels_[i]])
        
    #ax1.scatter(X[:,0], kmeans.labels_, label='k='+str(k))
    ax1.scatter(max_semimaj, max_mass,color=colorlist)

    ax1.set_yscale('log')
    ax1.set_xscale('log')
    ax1.set_xlabel('Maximum Semimajor Axis (MJup)')
    ax1.set_ylabel('Maximum Mass (AU)')

    plt.show()
#    ax1.clear()



    
