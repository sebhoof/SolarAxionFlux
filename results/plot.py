import numpy as np
import matplotlib.pyplot as plt

plt.rc('text', usetex=True)
plt.rc('text.latex', preamble=r'\usepackage{amsmath}\usepackage{amssymb}\usepackage{siunitx}\usepackage[cm]{sfmath}')

res1 = np.genfromtxt("primakoff.dat")
res2 = np.genfromtxt("compton.dat")
ref1 = np.genfromtxt("2013_redondo_primakoff.dat")
ref2 = np.genfromtxt("2013_redondo_compton.dat")

plt.plot(res1[:,0], 1.0e-4*50.0*res1[:,1]/1.0e20, 'r--')
plt.plot(res2[:,0], res2[:,1]/1.0e20, 'b--')
plt.plot(ref1[:,0], ref1[:,1]/1.0e20, 'ro', label=r'Primakoff')
plt.plot(ref2[:,0], 0.5*ref2[:,1]/1.0e20, 'bo', label=r'Compton')

plt.xlabel(r'Energy $\omega$ [keV]')
plt.ylabel(r'Axion flux $\mathrm{d}\Phi_a/\mathrm{d}\omega$ [\SI{e20}{\per\second\per\cm\squared\per\keV}]')

plt.legend()

plt.savefig("validation.pdf")
plt.show()
