import numpy as np
import matplotlib.pyplot as plt
from scipy import interpolate

plt.rc('text', usetex=True)
plt.rc('text.latex', preamble=r'\usepackage{amsmath}\usepackage{amssymb}\usepackage{siunitx}\usepackage[cm]{sfmath}\DeclareSIUnit\year{yr}')
conversion = 365.0*24.0*60.0*60.0*1.0e4*1.0e-20
res1 = np.genfromtxt("primakoff.dat")
res2 = np.genfromtxt("compton.dat")
res3 = np.genfromtxt("all_ff.dat")
res4 = np.genfromtxt("all_gaee.dat")
res5 = np.genfromtxt("metals.dat")
#corr = np.genfromtxt("weighted_compton.dat")
#weighted_compton = interpolate.interp1d(corr[:,0], corr[:,1], bounds_error=False, fill_value=0)
ref1 = np.genfromtxt("2013_redondo_primakoff.dat")
ref2 = np.genfromtxt("2013_redondo_compton.dat")
compton = interpolate.interp1d(ref2[:,0], ref2[:,1], bounds_error=False, fill_value=0)
ref3 = np.genfromtxt("2013_redondo_FF.dat")
ref4 = np.genfromtxt("2013_redondo_all.dat")

conv_fac = 1.0e-4/(365.0*24.0*60.0*60.0*1.0e10)

## Validation plots for axion-photon interactions

# Primakoff approximation [hep-ex/0702006] based on [astro-ph/0402114]
omega = np.linspace(0,10,300)
plt.plot(omega, 6.02*omega**2.481*np.exp(-omega/1.205),'b:', label=r'Primakoff approx. (BP04)')
plt.plot(ref1[:,0], conv_fac*(1.0e4/50.0)*ref1[:,1], 'r-', label=r'Primakoff (Redondo)')
plt.plot(res1[:,0], res1[:,1]/1.0e10, 'k--', label=r'Primakoff (AGSS09)')

plt.title(r'Axion-photon interactions, $g_{a\gamma\gamma} = \SI{e-10}{\GeV^{-1}}$, OP opacities')
plt.xlabel(r'Energy $\omega$ [keV]')
plt.ylabel(r'Axion flux $\mathrm{d}\Phi_a/\mathrm{d}\omega$ [\SI{e10}{\per\cm\squared\per\keV\per\s}]')
plt.xlim([0,10])
#plt.ylim([0,8])

plt.legend()

plt.savefig("validation_gagg.pdf")
plt.show()

## Validation plots for axion-electron interactions

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

plt.legend()

plt.savefig("validation_gaee.pdf")
plt.show()


gagg_bp98 = np.genfromtxt("gagg_bp98.dat")
gagg_b16gs98 = np.genfromtxt("gagg_b16gs98.dat")
gagg_b16agss09 = np.genfromtxt("gagg_b16agss09.dat")
gagg_agss09 = np.genfromtxt("gagg_agss09.dat")
gagg_agss09ph = np.genfromtxt("gagg_agss09ph.dat")
gagg_gs98 = np.genfromtxt("gagg_gs98.dat")
gagg_bp00 = np.genfromtxt("gagg_bp00.dat")
gagg_bp04 = np.genfromtxt("gagg_bp04.dat")
gagg_bs05op = np.genfromtxt("gagg_bs05op.dat")
gagg_bs05agsop = np.genfromtxt("gagg_bs05agsop.dat")

plt.plot(gagg_b16agss09[:,0], gagg_b16agss09[:,1]/1.0e10, 'r-', label=r'B16-AGSS09')
plt.plot(gagg_b16gs98[:,0], gagg_b16gs98[:,1]/1.0e10, '-', color='deepskyblue', label=r'B16-GS98')
plt.plot(gagg_agss09[:,0], gagg_agss09[:,1]/1.0e10, 'b:', label=r'AGSS09')
plt.plot(gagg_agss09ph[:,0], gagg_agss09ph[:,1]/1.0e10, 'g--', label=r'AGSS09ph')
plt.plot(gagg_gs98[:,0], gagg_gs98[:,1]/1.0e10, 'm-.', label=r'GS98')
plt.plot(gagg_bs05agsop[:,0], gagg_bs05agsop[:,1]/1.0e10, '-.', color='yellowgreen', label=r'BS05-AGSOP')
plt.plot(gagg_bs05op[:,0], gagg_bs05op[:,1]/1.0e10, '--', color='saddlebrown', label=r'BS05-OP')
plt.plot(gagg_bp04[:,0], gagg_bp04[:,1]/1.0e10, ':', color='mediumvioletred', label=r'BP04')
plt.plot(gagg_bp00[:,0], gagg_bp00[:,1]/1.0e10, '-', color='orange', label=r'BP00')
plt.plot(gagg_bp98[:,0], gagg_bp98[:,1]/1.0e10, '-.', color='gold', label=r'BP98')

plt.title(r'Axion-photon interactions, $g_{a\gamma\gamma} = \SI{e-10}{\GeV^{-1}}$, OP opacities')
plt.xlabel(r'Energy $\omega$ [keV]')
plt.ylabel(r'Axion flux $\mathrm{d}\Phi_a/\mathrm{d}\omega$ [\SI{e10}{\per\cm\squared\per\keV\per\s}]')
plt.xlim([0,10])
#plt.ylim([0,8])

plt.legend()

plt.savefig("solar_models_comp_gagg.pdf")
plt.show()

plt.plot(gagg_b16gs98[:,0], gagg_b16gs98[:,1]/gagg_b16agss09[:,1] - 1.0, '-', color='deepskyblue', label=r'B16-GS98')
plt.plot(gagg_agss09[:,0], gagg_agss09[:,1]/gagg_b16agss09[:,1] - 1.0, 'b:', label=r'AGSS09')
plt.plot(gagg_agss09ph[:,0], gagg_agss09ph[:,1]/gagg_b16agss09[:,1] - 1.0, 'g--', label=r'AGSS09ph')
plt.plot(gagg_gs98[:,0], gagg_gs98[:,1]/gagg_b16agss09[:,1] - 1.0, 'm-.', label=r'GS98')
plt.plot(gagg_bs05agsop[:,0], gagg_bs05agsop[:,1]/gagg_b16agss09[:,1] - 1.0, '-.', color='yellowgreen', label=r'BS05-AGSOP')
plt.plot(gagg_bs05op[:,0], gagg_bs05op[:,1]/gagg_b16agss09[:,1] - 1.0, '--', color='saddlebrown', label=r'BS05-OP')
plt.plot(gagg_bp04[:,0], gagg_bp04[:,1]/gagg_b16agss09[:,1] - 1.0, ':', color='mediumvioletred', label=r'BP04')
plt.plot(gagg_bp00[:,0], gagg_bp00[:,1]/gagg_b16agss09[:,1] - 1.0, '-', color='orange', label=r'BP00')
plt.plot(gagg_bp98[:,0], gagg_bp98[:,1]/gagg_b16agss09[:,1] - 1.0, '-.', color='gold', label=r'BP98')

plt.title(r'Axion-photon interactions, $g_{a\gamma\gamma}$, OP opacities')
plt.xlabel(r'Energy $\omega$ [keV]')
plt.ylabel(r'Relative deviation w.r.t. B16-AGSS09')
plt.xlim([0,10])
#plt.ylim([0,8])

plt.legend()

plt.savefig("solar_models_comp_gagg_relative.pdf")
plt.show()


gaee_bp98 = np.genfromtxt("gaee_bp98.dat")
gaee_b16gs98 = np.genfromtxt("gaee_b16gs98.dat")
gaee_b16agss09 = np.genfromtxt("gaee_b16agss09.dat")
gaee_agss09 = np.genfromtxt("gaee_agss09.dat")
gaee_agss09ph = np.genfromtxt("gaee_agss09ph.dat")
gaee_gs98 = np.genfromtxt("gaee_gs98.dat")
gaee_bp00 = np.genfromtxt("gaee_bp00.dat")
gaee_bp04 = np.genfromtxt("gaee_bp04.dat")
gaee_bs05op = np.genfromtxt("gaee_bs05op.dat")
gaee_bs05agsop = np.genfromtxt("gaee_bs05agsop.dat")

plt.plot(gaee_b16agss09[:,0], gaee_b16agss09[:,1]/1.0e8, 'r-', label=r'B16-AGSS09')
plt.plot(gaee_b16gs98[:,0], gaee_b16gs98[:,1]/1.0e8, '-', color='deepskyblue', label=r'B16-GS98')
plt.plot(gaee_agss09[:,0], gaee_agss09[:,1]/1.0e8, 'b:', label=r'AGSS09')
plt.plot(gaee_agss09ph[:,0], gaee_agss09ph[:,1]/1.0e8, 'g--', label=r'AGSS09ph')
plt.plot(gaee_gs98[:,0], gaee_gs98[:,1]/1.0e8, 'm-.', label=r'GS98')
plt.plot(gaee_bs05agsop[:,0], gaee_bs05agsop[:,1]/1.0e8, '-.', color='yellowgreen', label=r'BS05-AGSOP')
plt.plot(gaee_bs05op[:,0], gaee_bs05op[:,1]/1.0e8, '--', color='saddlebrown', label=r'BS05-OP')
plt.plot(gaee_bp04[:,0], gaee_bp04[:,1]/1.0e8, ':', color='mediumvioletred', label=r'BP04')
plt.plot(gaee_bp00[:,0], gaee_bp00[:,1]/1.0e8, '-', color='orange', label=r'BP00')
plt.plot(gaee_bp98[:,0], gaee_bp98[:,1]/1.0e8, '-.', color='gold', label=r'BP98')

plt.title(r'Axion-electron interactions, $g_{aee} = \num{e-13}$, OP opacities')
plt.xlabel(r'Energy $\omega$ [keV]')
plt.ylabel(r'Axion flux $\mathrm{d}\Phi_a/\mathrm{d}\omega$ [\SI{e8}{\per\cm\squared\per\keV\per\s}]')
plt.xlim([0,10])
plt.ylim([0,12])

plt.legend()

plt.savefig("solar_models_comp_gaee.pdf")
plt.show()

plt.plot(gaee_b16gs98[:,0], gaee_b16gs98[:,1]/gaee_b16agss09[:,1] - 1.0, '-', color='deepskyblue', label=r'B16-GS98')
plt.plot(gaee_agss09[:,0], gaee_agss09[:,1]/gaee_b16agss09[:,1] - 1.0, 'b:', label=r'AGSS09')
plt.plot(gaee_agss09ph[:,0], gaee_agss09ph[:,1]/gaee_b16agss09[:,1] - 1.0, 'g--', label=r'AGSS09ph')
plt.plot(gaee_gs98[:,0], gaee_gs98[:,1]/gaee_b16agss09[:,1] - 1.0, 'm-.', label=r'GS98')
plt.plot(gaee_bs05agsop[:,0], gaee_bs05agsop[:,1]/gaee_b16agss09[:,1] - 1.0, '-.', color='yellowgreen', label=r'BS05-AGSOP')
plt.plot(gaee_bs05op[:,0], gaee_bs05op[:,1]/gaee_b16agss09[:,1] - 1.0, '--', color='saddlebrown', label=r'BS05-OP')
plt.plot(gaee_bp04[:,0], gaee_bp04[:,1]/gaee_b16agss09[:,1] - 1.0, ':', color='mediumvioletred', label=r'BP04')
plt.plot(gaee_bp00[:,0], gaee_bp00[:,1]/gaee_b16agss09[:,1] - 1.0, '-', color='orange', label=r'BP00')
plt.plot(gaee_bp98[:,0], gaee_bp98[:,1]/gaee_b16agss09[:,1] - 1.0, '-.', color='gold', label=r'BP98')

plt.title(r'Axion-electron interactions, $g_{aee}$, OP opacities')
plt.xlabel(r'Energy $\omega$ [keV]')
plt.ylabel(r'Relative deviation w.r.t. B16-AGSS09')
plt.xlim([0,10])
#plt.ylim([0,10])

plt.legend()

plt.savefig("solar_models_comp_gaee_relative.pdf")
plt.show()


gaee_op = np.genfromtxt("gaee_OP.dat")
gaee_opas = np.genfromtxt("gaee_OPAS.dat")
gaee_ledcop = np.genfromtxt("gaee_LEDCOP.dat")
gaee_atomic = np.genfromtxt("gaee_ATOMIC.dat")

plt.plot(gaee_op[:,0], gaee_op[:,1]/1.0e8, 'r-', label=r'OP')
plt.plot(gaee_opas[gaee_opas[:,0]>=2.0,0], gaee_opas[gaee_opas[:,0]>=2.0,1]/1.0e8, 'b:', label=r'OPAS')
plt.plot(gaee_ledcop[:,0], gaee_ledcop[:,1]/1.0e8, 'g--', label=r'LEDCOP')
plt.plot(gaee_atomic[:,0], gaee_atomic[:,1]/1.0e8, 'm-.', label=r'ATOMIC')

plt.title(r'Axion-electron interactions, $g_{aee} = \num{e-13}$, AGSS09 model')
plt.xlabel(r'Energy $\omega$ [keV]')
plt.ylabel(r'Axion flux $\mathrm{d}\Phi_a/\mathrm{d}\omega$ [\SI{e8}{\per\cm\squared\per\keV\per\s}]')
plt.xlim([0,10])
plt.ylim([0,12])

plt.legend()

plt.savefig("opacity_codes_comp_gaee.pdf")
plt.show()


plt.plot(gaee_opas[gaee_opas[:,0]>=2.0,0], gaee_opas[gaee_opas[:,0]>=2.0,1]/gaee_op[gaee_opas[:,0]>=2.0,1] - 1.0, 'b:', label=r'OPAS')
plt.plot(gaee_ledcop[:,0], gaee_ledcop[:,1]/gaee_op[:,1] - 1.0, 'g--', label=r'LEDCOP')
plt.plot(gaee_atomic[:,0], gaee_atomic[:,1]/gaee_op[:,1] - 1.0, 'm-.', label=r'ATOMIC')

plt.title(r'Axion-electron interactions, $g_{aee}$, AGSS09 model')
plt.xlabel(r'Energy $\omega$ [keV]')
plt.ylabel(r'Relative deviation w.r.t. OP')
plt.xlim([0,10])
#plt.ylim([0,10])

plt.legend()

plt.savefig("opacity_codes_comp_gaee_relative.pdf")
plt.show()


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
# plt.plot(ref1[:,0], ref1[:,1]/1.0e20, 'r-', label=r'Primakoff')
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
# plt.show()
