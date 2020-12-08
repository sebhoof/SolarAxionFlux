import numpy as np
import matplotlib.pyplot as plt
from scipy import interpolate


plt.rc('text', usetex=True)
plt.rc('text.latex', preamble=r'\usepackage{amsmath}\usepackage{amssymb}\usepackage{siunitx}')

# Colours
col_b16agss09 = '#A50026'
col_b16gs98 = '#D73027'
col_agss09 = '#F46D43'
col_agss09ph = '#FDAE61'
col_ags05 = '#fEE090'
col_bs05agsop = '#FFFFBF'
col_bs05op = '#E0F3F8'
col_bp04 = '#ABD9E9'
col_bp00 = '#74ADD1'
col_bp98 = '#4575B4'
col_gs98 = '#313695'

def plot_setup(size=6,ratio=0.618):
    fig.set_size_inches(size,ratio*size)
    ax.tick_params(which='both', direction='in', bottom=True, top=True, left=True, right=True)
    ax.tick_params(which='major', length=6)
    ax.tick_params(which='minor', length=4)
    #plt.minorticks_on()

conversion = 365.0*24.0*60.0*60.0*1.0e4*1.0e-20
res1 = np.genfromtxt("primakoff.dat")
res2 = np.genfromtxt("compton.dat")
res3 = np.genfromtxt("all_ff.dat")
res4 = np.genfromtxt("all_gaee.dat")
res5 = np.genfromtxt("metals.dat")
res6 = np.genfromtxt("TP.dat")
res7 = np.genfromtxt("LP.dat")

#corr = np.genfromtxt("weighted_compton.dat")
#weighted_compton = interpolate.interp1d(corr[:,0], corr[:,1], bounds_error=False, fill_value=0)
common_path = "../data/benchmarks/"
ref1 = np.genfromtxt(common_path+"2013_redondo_primakoff.dat")
ref2 = np.genfromtxt(common_path+"2013_redondo_compton.dat")
compton = interpolate.interp1d(ref2[:,0], ref2[:,1], bounds_error=False, fill_value=0)
ref3 = np.genfromtxt(common_path+"2013_redondo_ff.dat")
ref4 = np.genfromtxt(common_path+"2013_redondo_all.dat")
ref5 = np.genfromtxt(common_path+"2020_giannotti_TP.dat")
ref6 = np.genfromtxt(common_path+"2020_giannotti_LP.dat")
ref7 = np.genfromtxt(common_path+"2020-o'hare.dat")
ref8 = np.genfromtxt(common_path+"2020_caputo_LP.dat")

conv_fac = 1.0e-4/(365.0*24.0*60.0*60.0*1.0e10)

## Validation plots for axion-photon interactions

# Primakoff approximation [hep-ex/0702006] based on [astro-ph/0402114]
omega = np.linspace(0,10,300)

fig, ax = plt.subplots()
plot_setup()
plt.plot(omega, 6.02*omega**2.481*np.exp(-omega/1.205),':', color=col_agss09, label=r'Primakoff approx. (BP04)')
plt.plot(ref1[:,0], conv_fac*(1.0e4/50.0)*ref1[:,1], '-', color=col_b16agss09, label=r'Primakoff (Redondo)')
plt.plot(res1[:,0], res1[:,1]/1.0e10, 'k--', label=r'Primakoff (AGSS09)')
plt.plot(res6[:,0], res6[:,1]/1.0e10, 'k--', label=r'TP (AGSS09)')

plt.title(r'Axion-photon interactions, $g_{a\gamma\gamma} = \SI{e-10}{\GeV^{-1}}$, OP opacities')
plt.xlabel(r'Energy $\omega$ [keV]')
plt.ylabel(r'Axion flux $\mathrm{d}\Phi_a/\mathrm{d}\omega$ [\SI{e10}{\per\cm\squared\per\keV\per\s}]')
plt.xlim([0,10])
#plt.ylim([0,8])

plt.legend(frameon=False)

plt.savefig("validation_gagg.pdf", bbox_inches='tight')
#plt.show()
plt.close()


fig, ax = plt.subplots()
plot_setup()
plt.plot(omega, 6.02*omega**2.481*np.exp(-omega/1.205),':', color=col_agss09, label=r'Primakoff approx. (BP04)')
plt.plot(ref1[:,0], conv_fac*(1.0e4/50.0)*ref1[:,1], '-', color=col_b16agss09, label=r'Primakoff (Redondo)')
plt.plot(res1[:,0], res1[:,1]/1.0e10, 'k--', label=r'Primakoff (AGSS09)')
plt.plot(res6[:,0], res6[:,1]/1.0e10, 'k-', label=r'TP (AGSS09)')
plt.plot(ref5[:,0], ref5[:,1]*4.0*1.4995, '-', color='green', label=r'TP (Giannotti)')#correct B conversion in giannotti result and adjust coupling constant

plt.title(r'Axion-photon interactions, $g_{a\gamma\gamma} = \SI{e-10}{\GeV^{-1}}$, OP opacities')
plt.xlabel(r'Energy $\omega$ [keV]')
plt.ylabel(r'Axion flux $\mathrm{d}\Phi_a/\mathrm{d}\omega$ [\SI{e10}{\per\cm\squared\per\keV\per\s}]')
plt.xlim([0.1,10])
plt.yscale('log')
plt.xscale('log')
#plt.ylim([0,8])

plt.legend(frameon=False)

plt.savefig("validation_Tplasmon.pdf", bbox_inches='tight')
plt.show()
plt.close()


fig, ax = plt.subplots()
plot_setup()
plt.plot(omega, 6.02*omega**2.481*np.exp(-omega/1.205),':', color=col_agss09, label=r'Primakoff approx. (BP04)')
plt.plot(ref1[:,0], conv_fac*(1.0e4/50.0)*ref1[:,1], '-', color=col_b16agss09, label=r'Primakoff (Redondo)')
plt.plot(res1[:,0], res1[:,1]/1.0e10, 'k--', label=r'Primakoff (AGSS09)')
plt.plot(res7[:,0], res7[:,1]/1.0e10*2.0, 'k-', label=r'LP (AGSS09)')
plt.plot(ref6[:,0], ref6[:,1]*4.0, '--', color='green', label=r'LP (Giannotti)') # correct coupling
plt.plot(ref7[:,0], ref7[:,1]/1.0e10*4.0/1.8, '--', color='orange', label=r'LP (O´Hare)') # correct coupling and angular average
plt.plot(ref8[:,0], ref8[:,1]/1.0e10*(3.0/5.0)**2, '--', color='gold', label=r'LP (Caputo)') #correct  field values

plt.title(r'Axion-photon interactions, $g_{a\gamma\gamma} = \SI{e-10}{\GeV^{-1}}$, OP opacities')
plt.xlabel(r'Energy $\omega$ [keV]')
plt.ylabel(r'Axion flux $\mathrm{d}\Phi_a/\mathrm{d}\omega$ [\SI{e10}{\per\cm\squared\per\keV\per\s}]')
plt.xlim([0.001,0.4])
#plt.yscale('log')
#plt.xscale('log')
plt.ylim([0.0,37])

plt.legend(frameon=False)

plt.savefig("validation_Lplasmon.pdf", bbox_inches='tight')
plt.show()
plt.close()

fig, ax = plt.subplots()
## Validation plots for axion-electron interactions
plot_setup()
plt.plot(ref2[:,0], 100.0*conv_fac*(0.5*ref2[:,1]), 'b-', label=r'Compton (Redondo)')
plt.plot(ref3[:,0], 100.0*conv_fac*ref3[:,1], 'm-', label=r'FF (Redondo)')
plt.plot(ref4[:,0], 1.0e11*ref4[:,1]*(1.0e-13/0.511e-10)**2/(24.0*60.0*60.0) - 100.0*conv_fac*(0.5*compton(ref4[:,0])), 'g-', label=r'All')
plt.plot(res2[:,0], res2[:,1]/1.0e8, 'k--', label=r'Compton (B16-AGSS09)')
plt.plot(res3[:,0], res3[:,1]/1.0e8, 'k--', label=r'FF (B16-AGSS09)')
plt.plot(res4[:,0], res4[:,1]/1.0e8, 'k--', label=r'All (B16-AGSS09)')
plt.plot(res5[:,0], res5[:,1]/1.0e8, 'k--', label=r'Metals (B16-AGSS09)')

plt.title(r'Axion-electron interactions, $g_{aee} = \num{e-13}$, OP opacities')
plt.xlabel(r'Energy $\omega$ [keV]')
plt.ylabel(r'Axion flux $\mathrm{d}\Phi_a/\mathrm{d}\omega$ [\SI{e8}{\per\cm\squared\per\keV\per\s}]')
plt.xlim([0,10])
plt.ylim([0,12])

plt.legend(ncol=2, frameon=False)

plt.savefig("validation_gaee.pdf")
#plt.show()
plt.close()
