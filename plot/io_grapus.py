#
#
# Python module for grapus output data
# Helpful global variables and functions 
# for reading, writing data
#
#


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

# Log Data

nlogcol=7
logcols = range(nlogcol)
logkeys = ['istar', 'mstar', 'mdisc', 'q_disc', 'rout', 'rfrag', 'nembryo']
loglabels = ['istar', 'mstar', 'mdisc', 'q_disc', 'rout', 'rfrag', 'nembryo']

logcoldict = dict(zip(logkeys,logcols))
loglabeldict = dict(zip(logkeys,loglabels))

logfmt = '%i %.8e % .8e %.8e %.8e %.8e %i'


# Colour scheme for plotting system portraits

 brown = '#663300'  # Brown Dwarfs
 red = 'red'        # Giant Planets
 blue = '#0099ff'   # Terrestrial Planets
