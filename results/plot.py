import numpy as np
import matplotlib.pyplot as plt
from scipy import interpolate

plt.rc('text', usetex=True)
plt.rc('text.latex', preamble=r'\usepackage{amsmath}\usepackage{amssymb}\usepackage{siunitx}\usepackage[cm]{sfmath}\DeclareSIUnit\year{yr}')

res1 = np.genfromtxt("primakoff.dat")
res2 = np.genfromtxt("compton.dat")
res3 = np.genfromtxt("all_ff.dat")
res4 = np.genfromtxt("all_gaee.dat")
ref1 = np.genfromtxt("2013_redondo_primakoff.dat")
ref2 = np.genfromtxt("2013_redondo_compton.dat")
compton = interpolate.interp1d(ref2[:,0], ref2[:,1], bounds_error=False, fill_value=0)
ref3 = np.genfromtxt("2013_redondo_FF.dat")
ref4 = np.genfromtxt("2013_redondo_all.dat")

plt.plot(ref1[:,0], ref1[:,1]/1.0e20, 'r-', label=r'Primakoff')
plt.plot(ref2[:,0], 0.5*ref2[:,1]/1.0e20, 'b-', label=r'Compton')
plt.plot(ref3[:,0], ref3[:,1]/1.0e20, 'm-', label=r'FF')
plt.plot(ref4[:,0], 365.0*1.0e4*0.1*ref4[:,1]*(1.0e-13/0.511e-10)**2 - 0.5*compton(ref4[:,0])/1.0e20, 'g-', label=r'Full axion-electron')
plt.plot(res1[:,0], 1.0e-4*50.0*res1[:,1]/1.0e20, 'k--')
plt.plot(res2[:,0], res2[:,1]/1.0e20, 'k--')
plt.plot(res3[:,0], res3[:,1]/1.0e20, 'k--')
plt.plot(res4[:,0], res4[:,1]/1.0e20, 'ko')

plt.xlabel(r'Energy $\omega$ [keV]')
plt.ylabel(r'Axion flux $\mathrm{d}\Phi_a/\mathrm{d}\omega$ [\SI{e20}{\per\m\squared\per\keV\per\year}]')
plt.xlim([0,10])

plt.legend()

plt.savefig("validation.pdf")
plt.show()
