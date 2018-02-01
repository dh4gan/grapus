# Written 26/1/17 by dh4gan
# Takes outputs from multiple runs of the grapus model and stitches them together
# Primarily, this is ensuring that the 'istar' variable is corrected

import numpy as np
import filefinder as ff
from sys import exit
from io_grapus import ninitialcol,nlogcol,nfinalcol,logcoldict,finalcoldict, initialfmt,finalfmt,logfmt

# Find three (sorted) file lists containing all initial, final and log files to be built


initialfiles = ff.find_sorted_local_input_fileset('*.initial')
finalfiles = ff.find_sorted_local_input_fileset('*.final')
logfiles = ff.find_sorted_local_input_fileset('*.log')

ninitial = len(initialfiles)
nfinal = len(finalfiles)
nlog = len(logfiles)

if ninitial!=nfinal or nfinal!=nlog:
    print 'Unequal file counts: '
    print 'Initial files ',ninitial 
    print 'Final files ',nfinal
    print 'Log files ',nlog
    exit()

# Define the total file names

totalinitialfile = 'total.initial'
totalfinalfile = 'total.final'
totallogfile = 'total.log'


# Begin with the initial files, as they are easy

print 'Stitching together .initial files'

first = True

for initialfile in initialfiles:
    
    # Read file 
    
    print 'Reading initial file ', initialfile
    initialdata = np.genfromtxt(initialfile)
    ninitialrow = initialdata.size/ninitialcol  

    print 'There are ',ninitialrow, ' rows'

    initialdata.reshape(ninitialrow,ninitialcol)
    
    if(first):
        totalinitialdata = initialdata
        first = False
    else:
        totalinitialdata = np.concatenate((totalinitialdata, initialdata),axis=0)
    

# Write data to total file

np.savetxt(totalinitialfile,totalinitialdata,fmt=initialfmt)

print 'Initial files stitched and saved to ',totalinitialfile

print '-----------------------------------'

print 'Stitching together .log files'

first = True

istarmax = 0

for logfile in logfiles:
    
    # Read file 
    
    print 'Reading log file ',logfile
    logdata = np.genfromtxt(logfile)
    nlogrow = logdata.size/nlogcol  

    print 'There are ',nlogrow, ' rows'

    logdata.reshape(nlogrow,nlogcol)
     
    # If this is the first run, then istar will be correct
    if(first):
        totallogdata = logdata
        laststar = logdata[-1,logcoldict['istar']]
        istarmax = laststar
        first = False
        
    else:
        # Otherwise, istar needs to be modified to accept the last largest value    
        print laststar
        logdata[:,logcoldict['istar']] = logdata[:,logcoldict['istar']] + laststar
        laststar = logdata[-1,logcoldict['istar']]
        print laststar
        totallogdata = np.concatenate((totallogdata, logdata),axis=0)
        
        print totallogdata.shape
        
        
# Write data to total file

np.savetxt(totallogfile,totallogdata,fmt=logfmt)

print 'Log files stitched and saved to ',totallogfile


print '-----------------------------------'

print 'Stitching together .final files'

first = True

istarmax = 0

for finalfile in finalfiles:
    
    # Read file 
    
    print 'Reading final file ',finalfile
    finaldata = np.genfromtxt(finalfile)
    nfinalrow = finaldata.size/nfinalcol  

    print 'There are ',nfinalrow, ' rows'

    finaldata.reshape(nfinalrow,nfinalcol)
     
    # If this is the first run, then istar will be correct
    if(first):
        totalfinaldata = finaldata
        laststar = finaldata[-1,finalcoldict['istar']]
        istarmax = laststar
        first = False
        
    # Otherwise, istar needs to be modified to accept the last largest value
    else:
        
        finaldata[:,finalcoldict['istar']] = finaldata[:,finalcoldict['istar']] + laststar
        laststar = finaldata[-1,finalcoldict['istar']]
        totalfinaldata = np.concatenate((totalfinaldata, finaldata),axis=0)
        
        print totalfinaldata.shape
        
        
# Write data to total file

np.savetxt(totalfinalfile,totalfinaldata,fmt=finalfmt)

print 'Final files stitched and saved to ',totalfinalfile
    
    