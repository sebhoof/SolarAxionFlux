import numpy as np
import matplotlib.pyplot as plt
from scipy import interpolate


plt.rc('text', usetex=True)
plt.rc('text.latex', preamble=r'\usepackage{amsmath}\usepackage{amssymb}\usepackage{siunitx}\DeclareSIUnit\year{yr}')
plt.rc('font', **{'family':'serif','size':10})
plt.rc('axes', labelsize=10)
plt.rc('xtick', **{'labelsize':10, 'major.pad':5})
plt.rc('ytick', **{'labelsize':10, 'major.pad':5})
plt.rc('legend', **{'fontsize':8, 'title_fontsize':10})
plt.rc('figure', titlesize=12)

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
plt.plot(ref5[:,0], ref5[:,1]*4.0, '-', color='green', label=r'TP (Giannotti)')

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

gagg_bp98 = np.genfromtxt("gagg_bp98.dat")
gagg_b16gs98 = np.genfromtxt("gagg_b16gs98.dat")
gagg_b16agss09 = np.genfromtxt("gagg_b16agss09.dat")
gagg_agss09 = np.genfromtxt("gagg_agss09.dat")
gagg_agss09ph = np.genfromtxt("gagg_agss09ph.dat")
gagg_ags05 = np.genfromtxt("gagg_ags05.dat")
gagg_gs98 = np.genfromtxt("gagg_gs98.dat")
gagg_bp00 = np.genfromtxt("gagg_bp00.dat")
gagg_bp04 = np.genfromtxt("gagg_bp04.dat")
gagg_bs05op = np.genfromtxt("gagg_bs05op.dat")
gagg_bs05agsop = np.genfromtxt("gagg_bs05agsop.dat")

fig, ax = plt.subplots()
plot_setup()
plt.plot(gagg_b16agss09[:,0], gagg_b16agss09[:,1]/1.0e10, '-', color=col_b16agss09, label=r'B16-AGSS09')
plt.plot(gagg_b16gs98[:,0], gagg_b16gs98[:,1]/1.0e10, '-', color=col_b16gs98, label=r'B16-GS98')
plt.plot(gagg_agss09[:,0], gagg_agss09[:,1]/1.0e10, ':', color=col_agss09, label=r'AGSS09')
plt.plot(gagg_agss09ph[:,0], gagg_agss09ph[:,1]/1.0e10, '--', color=col_agss09ph, label=r'AGSS09ph')
plt.plot(gagg_ags05[:,0], gagg_ags05[:,1]/1.0e10, '--', color=col_ags05, label=r'AGS05')
plt.plot(gagg_gs98[:,0], gagg_gs98[:,1]/1.0e10, '-.', color=col_gs98, label=r'GS98')
plt.plot(gagg_bs05agsop[:,0], gagg_bs05agsop[:,1]/1.0e10, '-.', color=col_bs05agsop, label=r'BS05-AGSOP')
plt.plot(gagg_bs05op[:,0], gagg_bs05op[:,1]/1.0e10, '--', color=col_bs05op, label=r'BS05-OP')
plt.plot(gagg_bp04[:,0], gagg_bp04[:,1]/1.0e10, ':', color=col_bp04, label=r'BP04')
plt.plot(gagg_bp00[:,0], gagg_bp00[:,1]/1.0e10, '-', color=col_bp00, label=r'BP00')
plt.plot(gagg_bp98[:,0], gagg_bp98[:,1]/1.0e10, '-.', color=col_bp98, label=r'BP98')

plt.title(r'Axion-photon interactions, $g_{a\gamma\gamma} = \SI{e-10}{\GeV^{-1}}$, OP opacities')
plt.xlabel(r'Energy $\omega$ [keV]')
plt.ylabel(r'Axion flux $\mathrm{d}\Phi_a/\mathrm{d}\omega$ [\SI{e10}{\per\cm\squared\per\keV\per\s}]')
plt.xlim([0,10])
#plt.ylim([0,8])

plt.legend(ncol=2, frameon=False)

plt.savefig("solar_models_comp_gagg.pdf")
#plt.show()
plt.close()

fig, ax = plt.subplots()
plot_setup()
plt.plot(gagg_b16gs98[:,0], gagg_b16gs98[:,1]/gagg_b16agss09[:,1] - 1.0, '-', color=col_b16gs98, label=r'B16-GS98')
plt.plot(gagg_agss09[:,0], gagg_agss09[:,1]/gagg_b16agss09[:,1] - 1.0, ':', color=col_agss09, label=r'AGSS09')
plt.plot(gagg_agss09ph[:,0], gagg_agss09ph[:,1]/gagg_b16agss09[:,1] - 1.0, '--', color=col_agss09ph, label=r'AGSS09ph')
plt.plot(gagg_gs98[:,0], gagg_gs98[:,1]/gagg_b16agss09[:,1] - 1.0, '-.', color=col_gs98, label=r'GS98')
plt.plot(gagg_ags05[:,0], gagg_ags05[:,1]/gagg_b16agss09[:,1] - 1.0, '--', color=col_ags05, label=r'AGS05')
plt.plot(gagg_bs05agsop[:,0], gagg_bs05agsop[:,1]/gagg_b16agss09[:,1] - 1.0, '-.', color=col_bs05agsop, label=r'BS05-AGSOP')
plt.plot(gagg_bs05op[:,0], gagg_bs05op[:,1]/gagg_b16agss09[:,1] - 1.0, '--', color=col_bs05op, label=r'BS05-OP')
plt.plot(gagg_bp04[:,0], gagg_bp04[:,1]/gagg_b16agss09[:,1] - 1.0, ':', color=col_bp04, label=r'BP04')
plt.plot(gagg_bp00[:,0], gagg_bp00[:,1]/gagg_b16agss09[:,1] - 1.0, '-', color=col_bp00, label=r'BP00')
plt.plot(gagg_bp98[:,0], gagg_bp98[:,1]/gagg_b16agss09[:,1] - 1.0, '-.', color=col_bp98, label=r'BP98')

plt.title(r'Axion-photon interactions, $g_{a\gamma\gamma}$, OP opacities')
plt.xlabel(r'Energy $\omega$ [keV]')
plt.ylabel(r'Relative deviation w.r.t. B16-AGSS09')
plt.xlim([0,10])
#plt.ylim([0,8])

plt.legend()

plt.savefig("solar_models_comp_gagg_relative.pdf")
#plt.show()
plt.close()

fig, ax = plt.subplots()
plot_setup(3,1)
plt.plot(gagg_b16gs98[:,0], gagg_b16gs98[:,1]/gagg_b16agss09[:,1] - 1.0, '-', color=col_b16gs98, label=r'B16-GS98')
plt.plot(gagg_agss09[:,0], gagg_agss09[:,1]/gagg_b16agss09[:,1] - 1.0, ':', color=col_agss09, label=r'AGSS09')
plt.plot(gagg_agss09ph[:,0], gagg_agss09ph[:,1]/gagg_b16agss09[:,1] - 1.0, '--', color=col_agss09ph, label=r'AGSS09ph')
plt.plot(gagg_gs98[:,0], gagg_gs98[:,1]/gagg_b16agss09[:,1] - 1.0, '-.', color=col_gs98, label=r'GS98')
plt.plot(gagg_ags05[:,0], gagg_ags05[:,1]/gagg_b16agss09[:,1] - 1.0, '--', color=col_ags05, label=r'AGS05')

ax.minorticks_on()

plt.xlabel(r'Energy $\omega$ [keV]')
plt.ylabel(r'Relative flux deviation')
plt.xlim([0,10])
#plt.ylim([-0.03,0.13])
plt.ylim([-0.05,0.20])

plt.legend(frameon=False)
plt.tight_layout()
plt.savefig("solar_models_comp_gagg_relative_relevant.pdf")
#plt.show()
plt.close()

gaee_bp98 = np.genfromtxt("gaee_bp98.dat")
gaee_b16gs98 = np.genfromtxt("gaee_b16gs98.dat")
gaee_b16agss09 = np.genfromtxt("gaee_b16agss09.dat")
gaee_agss09 = np.genfromtxt("gaee_agss09.dat")
gaee_agss09ph = np.genfromtxt("gaee_agss09ph.dat")
gaee_ags05 = np.genfromtxt("gaee_ags05.dat")
gaee_gs98 = np.genfromtxt("gaee_gs98.dat")
gaee_bp00 = np.genfromtxt("gaee_bp00.dat")
gaee_bp04 = np.genfromtxt("gaee_bp04.dat")
gaee_bs05op = np.genfromtxt("gaee_bs05op.dat")
gaee_bs05agsop = np.genfromtxt("gaee_bs05agsop.dat")

fig, ax = plt.subplots()
plot_setup()
plt.plot(gaee_b16agss09[:,0], gaee_b16agss09[:,1]/1.0e8, '-', color=col_b16agss09, label=r'B16-AGSS09')
plt.plot(gaee_b16gs98[:,0], gaee_b16gs98[:,1]/1.0e8, '-', color=col_b16gs98, label=r'B16-GS98')
plt.plot(gaee_agss09[:,0], gaee_agss09[:,1]/1.0e8, ':', color=col_agss09, label=r'AGSS09')
plt.plot(gaee_agss09ph[:,0], gaee_agss09ph[:,1]/1.0e8, '--', color=col_agss09ph, label=r'AGSS09ph')
plt.plot(gaee_gs98[:,0], gaee_gs98[:,1]/1.0e8, '-.', color=col_gs98, label=r'GS98')
plt.plot(gaee_ags05[:,0], gaee_ags05[:,1]/1.0e8, '--', color=col_ags05, label=r'AGS05')
plt.plot(gaee_bs05agsop[:,0], gaee_bs05agsop[:,1]/1.0e8, '-.', color=col_bs05agsop, label=r'BS05-AGSOP')
plt.plot(gaee_bs05op[:,0], gaee_bs05op[:,1]/1.0e8, '--', color=col_bs05op, label=r'BS05-OP')
plt.plot(gaee_bp04[:,0], gaee_bp04[:,1]/1.0e8, ':', color=col_bp04, label=r'BP04')
plt.plot(gaee_bp00[:,0], gaee_bp00[:,1]/1.0e8, '-', color=col_bp00, label=r'BP00')
plt.plot(gaee_bp98[:,0], gaee_bp98[:,1]/1.0e8, '-.', color=col_bp98, label=r'BP98')

plt.title(r'Axion-electron interactions, $g_{aee} = \num{e-13}$, OP opacities')
plt.xlabel(r'Energy $\omega$ [keV]')
plt.ylabel(r'Axion flux $\mathrm{d}\Phi_a/\mathrm{d}\omega$ [\SI{e8}{\per\cm\squared\per\keV\per\s}]')
plt.xlim([0,10])
plt.ylim([0,12])

plt.legend()

plt.savefig("solar_models_comp_gaee.pdf")
#plt.show()
plt.close()

fig, ax = plt.subplots()
plot_setup()
plt.plot(gaee_b16gs98[:,0], gaee_b16gs98[:,1]/gaee_b16agss09[:,1] - 1.0, '-', color=col_b16gs98, label=r'B16-GS98')
plt.plot(gaee_agss09[:,0], gaee_agss09[:,1]/gaee_b16agss09[:,1] - 1.0, ':', color=col_agss09, label=r'AGSS09')
plt.plot(gaee_agss09ph[:,0], gaee_agss09ph[:,1]/gaee_b16agss09[:,1] - 1.0, '--', color=col_agss09ph, label=r'AGSS09ph')
plt.plot(gaee_ags05[:,0], gaee_ags05[:,1]/gaee_b16agss09[:,1] - 1.0, '--', color=col_ags05, label=r'AGS05')
plt.plot(gaee_gs98[:,0], gaee_gs98[:,1]/gaee_b16agss09[:,1] - 1.0, '-.', color=col_gs98, label=r'GS98')
plt.plot(gaee_bs05agsop[:,0], gaee_bs05agsop[:,1]/gaee_b16agss09[:,1] - 1.0, '-.', color=col_bs05agsop, label=r'BS05-AGSOP')
plt.plot(gaee_bs05op[:,0], gaee_bs05op[:,1]/gaee_b16agss09[:,1] - 1.0, '--', color=col_bs05op, label=r'BS05-OP')
plt.plot(gaee_bp04[:,0], gaee_bp04[:,1]/gaee_b16agss09[:,1] - 1.0, ':', color=col_bp04, label=r'BP04')
plt.plot(gaee_bp00[:,0], gaee_bp00[:,1]/gaee_b16agss09[:,1] - 1.0, '-', color=col_bp00, label=r'BP00')
plt.plot(gaee_bp98[:,0], gaee_bp98[:,1]/gaee_b16agss09[:,1] - 1.0, '-.', color=col_bp98, label=r'BP98')

plt.title(r'Axion-electron interactions, $g_{aee}$, OP opacities')
plt.xlabel(r'Energy $\omega$ [keV]')
plt.ylabel(r'Relative deviation w.r.t. B16-AGSS09')
plt.xlim([0,10])
#plt.ylim([0,10])

plt.legend()

plt.savefig("solar_models_comp_gaee_relative.pdf")
#plt.show()
plt.close()

fig, ax = plt.subplots()
plot_setup(3, 1)
plt.plot(gaee_b16gs98[:,0], gaee_b16gs98[:,1]/gaee_b16agss09[:,1] - 1.0, '-', color=col_b16gs98, label=r'B16-GS98')
plt.plot(gaee_agss09[:,0], gaee_agss09[:,1]/gaee_b16agss09[:,1] - 1.0, ':', color=col_agss09, label=r'AGSS09')
plt.plot(gaee_agss09ph[:,0], gaee_agss09ph[:,1]/gaee_b16agss09[:,1] - 1.0, '--', color=col_agss09ph, label=r'AGSS09ph')
plt.plot(gaee_ags05[:,0], gaee_ags05[:,1]/gaee_b16agss09[:,1] - 1.0, '--', color=col_ags05, label=r'AGS05')
plt.plot(gaee_gs98[:,0], gaee_gs98[:,1]/gaee_b16agss09[:,1] - 1.0, '-.', color=col_gs98, label=r'GS98')

ax.minorticks_on()

plt.xlabel(r'Energy $\omega$ [keV]')
plt.ylabel(r'Relative flux deviation')
plt.xlim([0,10])
plt.ylim([-0.05,0.20])
plt.legend(frameon=False)

plt.tight_layout()
plt.savefig("solar_models_comp_gaee_relative_relevant.pdf")
#plt.show()
plt.close()


fig, ax = plt.subplots()
plot_setup(5)

plt.plot(gaee_b16agss09[:,0], gaee_b16agss09[:,1]/1.0e8, '-', color=col_b16agss09, label=r'B16-AGSS09')
plt.plot(gaee_b16gs98[:,0], gaee_b16gs98[:,1]/1.0e8, '-', color=col_b16gs98, label=r'B16-GS98')
plt.plot(gaee_agss09[:,0], gaee_agss09[:,1]/1.0e8, ':', color=col_agss09, label=r'AGSS09')
plt.plot(gaee_agss09ph[:,0], gaee_agss09ph[:,1]/1.0e8, '--', color=col_agss09ph, label=r'AGSS09ph')
plt.plot(gaee_ags05[:,0], gaee_ags05[:,1]/1.0e8, '--', color=col_ags05, label=r'AGS05')
plt.plot(gaee_bs05agsop[:,0], gaee_bs05agsop[:,1]/1.0e8, '-.', color=col_bs05agsop, label=r'BS05-AGSOP')
plt.plot(gaee_bs05op[:,0], gaee_bs05op[:,1]/1.0e8, '--', color=col_bs05op, label=r'BS05-OP')
plt.plot(gaee_bp04[:,0], gaee_bp04[:,1]/1.0e8, ':', color=col_bp04, label=r'BP04')
plt.plot(gaee_bp00[:,0], gaee_bp00[:,1]/1.0e8, '-', color=col_bp00, label=r'BP00')
plt.plot(gaee_bp98[:,0], gaee_bp98[:,1]/1.0e8, '-.', color=col_bp98, label=r'BP98')
plt.plot(gaee_gs98[:,0], gaee_gs98[:,1]/1.0e8, '-.', color=col_gs98, label=r'GS98')

plt.plot(gagg_b16agss09[:,0], gagg_b16agss09[:,1]/1.0e10, '-', color=col_b16agss09)
plt.plot(gagg_b16gs98[:,0], gagg_b16gs98[:,1]/1.0e10, '-', color=col_b16gs98)
plt.plot(gagg_agss09[:,0], gagg_agss09[:,1]/1.0e10, ':', color=col_agss09)
plt.plot(gagg_agss09ph[:,0], gagg_agss09ph[:,1]/1.0e10, '--', color=col_agss09ph)
plt.plot(gagg_ags05[:,0], gagg_ags05[:,1]/1.0e10, '--', color=col_ags05)
plt.plot(gagg_bs05agsop[:,0], gagg_bs05agsop[:,1]/1.0e10, '-.', color=col_bs05agsop)
plt.plot(gagg_bs05op[:,0], gagg_bs05op[:,1]/1.0e10, '--', color=col_bs05op)
plt.plot(gagg_bp04[:,0], gagg_bp04[:,1]/1.0e10, ':', color=col_bp04)
plt.plot(gagg_bp00[:,0], gagg_bp00[:,1]/1.0e10, '-', color=col_bp00)
plt.plot(gagg_bp98[:,0], gagg_bp98[:,1]/1.0e10, '-.', color=col_bp98)
plt.plot(gagg_gs98[:,0], gagg_gs98[:,1]/1.0e10, '-.', color=col_gs98)

ax.minorticks_on()
ax.tick_params(axis='y', which='minor', left=False, right=False)
plt.xlabel(r'Energy $\omega$ [keV]')
plt.ylabel(r'Axion flux $\mathrm{d}\Phi_a/\mathrm{d}\omega$ [\SI{e10}{\per\cm\squared\per\keV\per\s}]')
plt.xlim([0,10])
#plt.ylim([0,8])

plt.legend(ncol=1, frameon=False, handlelength=3, handletextpad=0.8, labelspacing=0.5)

plt.savefig("solar_models_comp_both.pdf", bbox_inches='tight')
#plt.show()
plt.close()

gaee_op = np.genfromtxt("gaee_OP.dat")
gaee_opas = np.genfromtxt("gaee_OPAS.dat")
gaee_ledcop = np.genfromtxt("gaee_LEDCOP.dat")
gaee_atomic = np.genfromtxt("gaee_ATOMIC.dat")


fig, ax = plt.subplots()
plot_setup(3,1)
plt.plot(gaee_op[:,0], gaee_op[:,1]/1.0e8, '-', color=col_b16agss09, label=r'OP')
plt.plot(gaee_opas[gaee_opas[:,0]>=2.0,0], gaee_opas[gaee_opas[:,0]>=2.0,1]/1.0e8, ':', color=col_agss09, label=r'OPAS')
plt.plot(gaee_ledcop[:,0], gaee_ledcop[:,1]/1.0e8, '--', color=col_agss09ph, label=r'LEDCOP')
plt.plot(gaee_atomic[:,0], gaee_atomic[:,1]/1.0e8, '-.', color=col_gs98, label=r'ATOMIC')

#plt.title(r'Axion-electron interactions, $g_{aee} = \num{e-13}$, AGSS09 model')

ax.minorticks_on()
ax.tick_params(axis='y', which='minor', left=False, right=False)

plt.xlabel(r'Energy $\omega$ [keV]')
plt.ylabel(r'Axion flux $\mathrm{d}\Phi_a/\mathrm{d}\omega$ [\SI{e8}{\per\cm\squared\per\keV\per\s}]')
plt.xlim([0,10])
plt.ylim([0,12])

plt.legend(frameon=False)

plt.tight_layout()

plt.savefig("opacity_codes_comp_gaee.pdf")
#plt.show()
plt.close()


fig, ax = plt.subplots()
plot_setup(3,1)
plt.plot(gaee_opas[gaee_opas[:,0]>=2.0,0], gaee_opas[gaee_opas[:,0]>=2.0,1]/gaee_op[gaee_opas[:,0]>=2.0,1] - 1.0, '-', color=col_agss09, label=r'OPAS')
plt.plot(gaee_ledcop[:,0], gaee_ledcop[:,1]/gaee_op[:,1] - 1.0, '--', color=col_agss09ph, label=r'LEDCOP')
plt.plot(gaee_atomic[:,0], gaee_atomic[:,1]/gaee_op[:,1] - 1.0, '-.', color=col_gs98, label=r'ATOMIC')

#plt.title(r'Axion-electron interactions, $g_{aee}$, AGSS09 model')
plt.xlabel(r'Energy $\omega$ [keV]')
plt.ylabel(r'Relative flux deviations')
plt.xlim([0,10])
plt.ylim([-0.75,0.75])

ax.minorticks_on()

plt.legend(frameon=False)

plt.tight_layout()

plt.savefig("opacity_codes_comp_gaee_relative.pdf")
#plt.show()
plt.close()


# res1 = np.genfromtxt("primakoff.dat")
# res2 = np.genfromtxt("compton.dat")
# res3 = np.genfromtxt("all_ff.dat")
# res4 = np.genfromtxt("metals.dat")
#
# ref1 = np.genfromtxt("2013_redondo_primakoff.dat")
# ref2 = np.genfromtxt("2013_redondo_compton.dat")
# ref3 = np.genfromtxt("2013_redondo_FF.dat")
# ref4 = np.genfromtxt("2013_redondo_all.dat")
# #compton fit
# omega=np.linspace(0,10,1000)
# plt.plot(omega,2.0e-2*50*omega**2.45*np.exp(-0.829*omega))
#
# plt.plot(ref1[:,0], ref1[:,1]/1.0e20, '-', color=col_b16agss09, label=r'Primakoff')
# plt.plot(ref2[:,0], ref2[:,1]/1.0e20, 'b-', label=r'Compton')
# plt.plot(ref3[:,0], ref3[:,1]/1.0e20, 'm-', label=r'FF')
# plt.plot(ref4[:,0], 365.0*1.0e4*0.1*ref4[:,1]*(1.0e-13/0.511e-10)**2 , 'g-', label=r'Full axion-electron')
#
# plt.plot(res1[:,0], 1.0e-4*50.0*res1[:,1]/1.0e20, 'k--')
# plt.plot(res2[:,0], 2.0*res2[:,1]/1.0e20, 'k--')
# plt.plot(res3[:,0], res3[:,1]/1.0e20, 'k--')
# plt.plot(res4[:,0], (res4[:,1]+res3[:,1]+2.0*res2[:,1])/1.0e20, 'k--')
#
# plt.xlabel(r'Energy $\omega$ [keV]')
# plt.ylabel(r'Axion flux $\mathrm{d}\Phi_a/\mathrm{d}\omega$ [\SI{e20}{\per\m\squared\per\keV\per\year}]')
# plt.xlim([0,10])
#
# plt.legend()
#
# plt.savefig("copy_Javi.pdf")
# #plt.show()
