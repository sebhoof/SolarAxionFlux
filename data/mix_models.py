import numpy as np
import pylab as plt
import time
from scipy.interpolate import interp1d
import pexpect
from math import ceil
from scipy.integrate import dblquad



'''
data98=np.genfromtxt('model_gs98.dat',skip_header=0,skip_footer=2)
data09= np.genfromtxt('SolarModel_AGSS09ph.dat',skip_header=0,skip_footer=2)
print (len(data09[:,0]))
print (len(data98[:1967,0]))
print (np.shape(data09[:,:]))
newdata= np.empty((1967,35))
newdata[:,0:6]=data09[:,0:6]
newdata[:,6:35]=data98[:1967,6:35]
np.savetxt('model_Javi.txt',newdata)
'''

#sunRed= np.genfromtxt('/Users/lennertthormaehlen/Documents/Uni/Promotion/axionflux/redondo/OP_AR/OPAtables/solarmodelREDUCED_fine.dat')
solarmodel= np.genfromtxt('/Users/lennertthormaehlen/Documents/Uni/Promotion/axionflux/redondo/solardata/solarmodel.txt')
data98=np.genfromtxt('/Users/lennertthormaehlen/Documents/Uni/Promotion/axionflux/redondo/solardata/model_gs98.dat',skip_header=20,skip_footer=2)
solarmodelred= np.genfromtxt('/Users/lennertthormaehlen/Documents/Uni/Promotion/axionflux/redondo/OP_AR/OPAtables/solarmodelREDUCED_fine.dat')
print(solarmodelred)

print(np.shape(solarmodel))
print(np.shape(data98))

newdata = np.zeros((2000,35))
newdata[:,0:2] = data98[:,0:2]
newdata[:,5] = data98[:,5]
newdata[:,9:35] = data98[:,9:35]

X = interp1d(solarmodel[:,0],solarmodel[:,3])
newdata[:,6] = X(data98[:,1])
newdata[:,7] = 1.0-newdata[:,6]

temperature = interp1d(solarmodel[:,0],solarmodel[:,2])
newdata[:,2] = temperature(data98[:,1])/8.6173303e-5

density = interp1d(solarmodel[:,0],solarmodel[:,1])
newdata[:,3] = density(data98[:,1])*1.0e-9/(197.327053e-10**3.0)*1.660539*1e-24

ne = interp1d(solarmodel[:,0],solarmodel[:,4])
radiationpress = 4.0/3.0*5.678e-15*newdata[:,2]**4.0
newdata[:,4] = (ne(data98[:,1])+(3.0*newdata[:,6]+1.0)/4.0*density(data98[:,1]))*1.0e-9/(197.327053e-10**3.0)*newdata[:,2]*1.381e-16+radiationpress

newdatared = np.zeros((len(solarmodelred[:,0]),35))
newdatared[:,1] = solarmodelred[:,0]
#X and Y
newdatared[:,6] = solarmodelred[:,3]
newdatared[:,7] = 1.0 -solarmodelred[:,3]
#temperature
newdatared[:,2] = solarmodelred[:,2]/8.6173303e-5
#density
newdatared[:,3] = solarmodelred[:,1]*1.0e-9/(197.327053e-10**3.0)*1.660539*1e-24
#pressure
nered = interp1d(solarmodel[:,0],solarmodel[:,4])
radiationpressred = 4.0/3.0*5.678e-15*newdatared[:,2]**4.0
newdatared[:,4] = (ne(solarmodelred[:,0])+(3.0*newdatared[:,6]+1.0)/4.0*density(solarmodelred[:,0]))*1.0e-9/(197.327053e-10**3.0)*newdatared[:,2]*1.381e-16+radiationpressred

modelinterp = interp1d(newdata[:,1],newdata[:,9:35],axis=0,bounds_error=False, fill_value='extrapolate')
newdatared[:,9:35] = modelinterp(solarmodelred[:,0])
print(newdatared)
np.savetxt('model_Javi.txt',newdata)
np.savetxt('model_Javi_red.txt',newdatared)
