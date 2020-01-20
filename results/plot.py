import numpy as np
import matplotlib.pyplot as plt
from scipy import interpolate

plt.rc('text', usetex=True)
plt.rc('text.latex', preamble=r'\usepackage{amsmath}\usepackage{amssymb}\usepackage{siunitx}\usepackage[cm]{sfmath}\DeclareSIUnit\year{yr}')

res1 = np.genfromtxt("primakoff.dat")
res2 = np.genfromtxt("compton.dat")
res3 = np.genfromtxt("all_ff.dat")
res4 = np.genfromtxt("all_gaee.dat")
#corr = np.genfromtxt("weighted_compton.dat")
#weighted_compton = interpolate.interp1d(corr[:,0], corr[:,1], bounds_error=False, fill_value=0)
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
plt.plot(res4[:,0], res4[:,1]/1.0e20, 'k--')

plt.xlabel(r'Energy $\omega$ [keV]')
plt.ylabel(r'Axion flux $\mathrm{d}\Phi_a/\mathrm{d}\omega$ [\SI{e20}{\per\m\squared\per\keV\per\year}]')
plt.xlim([0,10])

plt.legend()

plt.savefig("validation.pdf")
plt.show()

int_ref1 = interpolate.interp1d(ref1[:,0], ref1[:,1], bounds_error=False, fill_value=np.nan)
int_ref2 = interpolate.interp1d(ref2[:,0], 0.5*ref2[:,1], bounds_error=False, fill_value=np.nan)
int_ref3 = interpolate.interp1d(ref3[:,0], ref3[:,1], bounds_error=False, fill_value=np.nan)
int_ref4 = interpolate.interp1d(ref4[:,0], 365.0*1.0e4*0.1*ref4[:,1]*(1.0e-13/0.511e-10)**2 - 0.5*compton(ref4[:,0])/1.0e20, bounds_error=False, fill_value=np.nan)
int_res1 = interpolate.interp1d(res1[:,0], 1.0e-4*50.0*res1[:,1], bounds_error=False, fill_value=np.nan)
int_res2 = interpolate.interp1d(res2[:,0], res2[:,1], bounds_error=False, fill_value=np.nan)
int_res3 = interpolate.interp1d(res3[:,0], res3[:,1], bounds_error=False, fill_value=np.nan)
int_res4 = interpolate.interp1d(res4[:,0], res4[:,1]/1.0e20, bounds_error=False, fill_value=np.nan)

ergs = res1[:,0]

plt.plot(ergs, int_res1(ergs)/int_ref1(ergs) - 1.0, 'r-', label=r'Primakoff')
plt.plot(ergs, int_res2(ergs)/int_ref2(ergs) - 1.0, 'b-', label=r'Compton')
plt.plot(ergs, int_res3(ergs)/int_ref3(ergs) - 1.0, 'm-', label=r'FF')
plt.plot(ergs, int_res4(ergs)/int_ref4(ergs) - 1.0, 'g-', label=r'Full axion-electron')

plt.xlabel(r'Energy $\omega$ [keV]')
plt.ylabel(r'Realtive deviation from reference')
plt.xlim([0,10])

plt.legend()

plt.savefig("discrepancies.pdf")
plt.show()
