from pysb import *
from pysb.simulator import *
import numpy as np
import matplotlib.pyplot as plt
import os

# test

Model()

Monomer('TGFB', ['r'])
Monomer('TGFBRII', ['l', 'smad', 'avb3'])
Monomer('gITGB3')
Monomer('mITGB3')
Monomer('aVb3', ['tgfbrii'])
Monomer('SMAD', ['tgfbrii', 'state'], {'state': ['U', 'P', 'SIS3']})
Monomer('RhoA', ['state'], {'state': ['I', 'A']})
Monomer('ROCK', ['state'], {'state': ['I', 'A']})
Monomer('MEK', ['state'], {'state': ['U', 'P']})
Monomer('p38MAPK', ['state'], {'state': ['U', 'P', 'SB202190']})
Monomer('gGli2')
Monomer('mGli2')
Monomer('Gli2', ['gant58'])
Monomer('gPTHrP')
Monomer('mPTHrP')
Monomer('PTHrP')

Initial(TGFB(r=None), Parameter('TGFB_init', 100))
Initial(TGFBRII(l=None, smad=None, avb3=None), Parameter('TGFBRII_init', 1000))
Initial(gITGB3(), Parameter('gITGB3_init', 1))
Initial(aVb3(tgfbrii=None), Parameter('aVb3_init', 10))
Initial(SMAD(tgfbrii=None, state='U'), Parameter('SMAD_init', 100))
Initial(RhoA(state='I'), Parameter('RhoA_init', 100))
Initial(ROCK(state='I'), Parameter('ROCK_init', 100))
Initial(MEK(state='U'), Parameter('MEK_init', 100))
Initial(p38MAPK(state='U'), Parameter('p38MAPK_init', 100))
Initial(gGli2(), Parameter('gGli2_init', 1))
Initial(gPTHrP(), Parameter('gPTHrP_init', 1))

# Observable('TGFB_tot', TGFB())
# Observable('TGFBRII_TGFB', TGFBRII(l=1) % TGFB(r=1))
# Observable('TGFBRII_aVb3', TGFBRII(avb3=1) % aVb3(tgfbrii=1))
Observable('SMAD_p', SMAD(tgfbrii=None, state='P'))
Observable('RhoA_a', RhoA(state='A'))
# # Observable('ROCK_a', ROCK(state='A'))
# # Observable('MEK_p', MEK(state='P'))
Observable('p38MAPK_p', p38MAPK(state='P'))
# Observable('aVb3_tot', aVb3())
# # Observable('Gli2_tot', Gli2())
Observable('Gli2_free', Gli2(gant58=None))
# Observable('PTHrP_tot', PTHrP())

# TGFB->SMAD pathway
Parameter('kf_tgfbrii_lig_bind', 10)
Parameter('kr_tgfbrii_lig_bind', 100)
Parameter('kf_tgfbrii_smad_bind', 1)
Parameter('kr_tgfbrii_smad_bind', 10)
Parameter('k_smad_phospho', 10)
Parameter('k_smad_dephospho', 10)
Rule('TGFBRII_lig_bind', TGFBRII(l=None) + TGFB(r=None) | TGFBRII(l=1) % TGFB(r=1), kf_tgfbrii_lig_bind,
     kr_tgfbrii_lig_bind)
Rule('TGFBRII_SMAD_bind', TGFBRII(l=ANY, smad=None) + SMAD(tgfbrii=None) >> TGFBRII(l=ANY, smad=1) % SMAD(tgfbrii=1),
     kf_tgfbrii_smad_bind)
Rule('TGFBRII_SMAD_unbind', TGFBRII(smad=1) % SMAD(tgfbrii=1) >> TGFBRII(smad=None) + SMAD(tgfbrii=None),
     kr_tgfbrii_smad_bind)
Rule('SMAD_phosphorylation', TGFBRII(l=ANY, smad=1) % SMAD(tgfbrii=1, state='U') >>
     TGFBRII(l=ANY, smad=1) % SMAD(tgfbrii=1, state='P'), k_smad_phospho)
Rule('SMAD_dephosphorylation', SMAD(tgfbrii=None, state='P') >> SMAD(tgfbrii=None, state='U'), k_smad_dephospho)

# TGFBRII-ITGB3->p38/MAPK pathway
Parameter('kf_tgfbrii_avb3_bind', 1)  # ***
Parameter('kr_tgfbrii_avb3_bind', 100)
Parameter('kact_basal_rhoa', 1e-3)
Parameter('kact_induced_rhoa', 10*kact_basal_rhoa.value)
Parameter('kdeact_rhoa', 10)
Parameter('kact_rock', 1e-2)
Parameter('kdeact_rock', 10)
Parameter('kact_mek', 1e-2)
Parameter('kdeact_mek', 10)
Parameter('kact_p38mapk', 1e-2)
Parameter('kdeact_p38mapk', 10)
Rule('TGFBRII_aVb3_bind', TGFBRII(avb3=None) + aVb3(tgfbrii=None) | TGFBRII(avb3=1) % aVb3(tgfbrii=1),
     kf_tgfbrii_avb3_bind, kr_tgfbrii_avb3_bind)
Rule('RhoA_basal_activation', aVb3(tgfbrii=None) + RhoA(state='I') >> aVb3(tgfbrii=None) + RhoA(state='A'),
     kact_basal_rhoa)
Rule('RhoA_induced_activation', aVb3(tgfbrii=ANY) + RhoA(state='I') >> aVb3(tgfbrii=ANY) + RhoA(state='A'),
     kact_induced_rhoa)
Rule('RhoA_deactivation', RhoA(state='A') >> RhoA(state='I'), kdeact_rhoa)
Rule('ROCK_activation', RhoA(state='A') + ROCK(state='I') >> RhoA(state='A') + ROCK(state='A'), kact_rock)
Rule('ROCK_deactivation', ROCK(state='A') >> ROCK(state='I'), kdeact_rock)
Rule('MEK_activation', ROCK(state='A') + MEK(state='U') >> ROCK(state='A') + MEK(state='P'), kact_mek)
Rule('MEK_deactivation', MEK(state='P') >> MEK(state='U'), kdeact_mek)
Rule('p38MAPK_activation', MEK(state='P') + p38MAPK(state='U') >> MEK(state='P') + p38MAPK(state='P'), kact_p38mapk)
Rule('p38MAPK_deactivation', p38MAPK(state='P') >> p38MAPK(state='U'), kdeact_p38mapk)

# Gli2 expression by SMAD + p38/MAPK
Parameter('ktr1_gli2', 100)
Parameter('ktr2_gli2', 1)
Expression('ktr_gli2_smad_p38mapk', ktr1_gli2 * SMAD_p**2 * p38MAPK_p**2 / (ktr2_gli2**4 + SMAD_p**2 * p38MAPK_p**2))
Parameter('ktl_gli2', 100)
Parameter('kdeg_mGli2', 1)
Parameter('kdeg_pGli2', 1)
Rule('Gli2_transcription', gGli2() >> gGli2() + mGli2(), ktr_gli2_smad_p38mapk)
Rule('Gli2_translation', mGli2() >> mGli2() + Gli2(gant58=None), ktl_gli2)
Rule('mGli2_degradation', mGli2() >> None, kdeg_mGli2)
Rule('pGli2_degradation', Gli2(gant58=None) >> None, kdeg_pGli2)

# ITGB3 expression by Gli2
Parameter('ktr_itgb3_basal', 0.1)
Parameter('ktr1_itgb3', 100)
Parameter('ktr2_itgb3', 10)
Expression('ktr_itgb3_gli2', ktr1_itgb3 * Gli2_free**2 / (ktr2_itgb3**2 + Gli2_free**2))
Parameter('ktl_itgb3', 100)
Parameter('kdeg_mITGB3', 1)
Parameter('kdeg_avb3', 1)
Rule('ITGB3_transcription_basal', gITGB3() >> gITGB3() + mITGB3(), ktr_itgb3_basal)
Rule('ITGB3_transcription', gITGB3() >> gITGB3() + mITGB3(), ktr_itgb3_gli2)
Rule('ITGB3_translation', mITGB3() >> mITGB3() + aVb3(tgfbrii=None), ktl_itgb3)
Rule('mITGB3_degradation', mITGB3() >> None, kdeg_mITGB3)
Rule('aVb3_degradation', aVb3(tgfbrii=None) >> None, kdeg_avb3)

# PTHrP expression by Gli2
Parameter('ktr1_pthrp', 100)
Parameter('ktr2_pthrp', 1000)
Expression('ktr_pthrp_gli2', ktr1_pthrp * Gli2_free**2 / (ktr2_pthrp**2 + Gli2_free**2))
Parameter('ktl_pthrp', 100)
Parameter('kdeg_mPTHrP', 1)
Parameter('kdeg_pPTHrP', 1)
Rule('PTHrP_transcription', gPTHrP() >> gPTHrP() + mPTHrP(), ktr_pthrp_gli2)
Rule('PTHrP_translation', mPTHrP() >> mPTHrP() + PTHrP(), ktl_pthrp)
Rule('mPTHrP_degradation', mPTHrP() >> None, kdeg_mPTHrP)
Rule('pPTHrP_degradation', PTHrP() >> None, kdeg_pPTHrP)

# TGFB flush
Parameter('k_tgfb_flush', 0)
Rule('TGFB_flush', TGFB(r=None) >> None, k_tgfb_flush)
########################

# 1D11 binds TGFb
Monomer('x1D11', ['tgfb'])
Initial(x1D11(tgfb=None), Parameter('x1D11_init', 0))
# Observable('x1D11_TGFB', TGFB(r=1) % x1D11(tgfb=1))
Parameter('kf_1d11_tgfb_bind', 1000)
Parameter('kr_1d11_tgfb_bind', 10)
Rule('x1D11_binds_TGFB', TGFB(r=None) + x1D11(tgfb=None) | TGFB(r=1) % x1D11(tgfb=1), kf_1d11_tgfb_bind,
     kr_1d11_tgfb_bind)
###########################

# GANT58 binds Gli2
Monomer('GANT58', ['gli2'])
Initial(GANT58(gli2=None), Parameter('GANT58_init', 0))
# Observable('GANT58_Gli2', GANT58(gli2=1) % Gli2(gant58=1))
Parameter('kf_gant58_gli2_bind', 1000)
Parameter('kr_gant58_gli2_bind', 10)
Rule('GANT58_binds_Gli2', GANT58(gli2=None) + Gli2(gant58=None) | GANT58(gli2=1) % Gli2(gant58=1), kf_gant58_gli2_bind,
     kr_gant58_gli2_bind)
#############################

# SMAD~U + SIS3 -> SMAD~SIS3 + SIS3
# SMAD~SIS3 -> SMAD~U
Monomer('SIS3')
Initial(SIS3(), Parameter('SIS3_init', 0))
Parameter('kf_smad_sis3', 1)
Parameter('kr_smad_sis3', 100)
Rule('SMAD_U_to_SMAD_SIS3', SMAD(state='U') + SIS3() >> SMAD(state='SIS3') + SIS3(), kf_smad_sis3)
Rule('SMAD_SIS3_to_SMAD_U', SMAD(state='SIS3') >> SMAD(state='U'), kr_smad_sis3)
#############################

# p38MAPK~U + SB202190 -> p38MAPK~SB202190 + SB202190
# SMAD~SB202190 -> SMAD~U
Monomer('SB202190')
Initial(SB202190(), Parameter('SB202190_init', 0))
Parameter('kf_p38mapk_sb202190', 1)
Parameter('kr_p38mapk_sb202190', 100)
Rule('p38MAPK_U_to_p38MAPK_SB202190', p38MAPK(state='U') + SB202190() >> p38MAPK(state='SB202190') + SB202190(),
     kf_p38mapk_sb202190)
Rule('p38MAPK_SB202190_to_p38MAPK_U', p38MAPK(state='SB202190') >> p38MAPK(state='U'), kr_p38mapk_sb202190)
#############################

# EQUILIBRATE WITH MECHANICAL FORCES
tspan1 = np.linspace(0, 40, 401)
sim = ScipyOdeSimulator(model, verbose=True)
x = sim.run(tspan=tspan1)  # , param_values = {'kf_tgfbrii_avb3_bind' : 0})

# BIFURCATION
# tspan_equil = np.linspace(0,1000,10001)
# KD = np.append(np.arange(1e2, 5e3+1, 100), np.arange(1e4, 1e6+1e4, 1e4))
# plt.figure()
# for color in ['b', 'r']:
#     if color == 'b':
#         x = sim.run(tspan = tspan_equil, param_values= {'kr_tgfbrii_avb3_bind' : 1e2})
#     else:
#         x = sim.run(tspan = tspan_equil, param_values= {'kf_tgfbrii_avb3_bind' : 0})
#     initials = x.species[-1]
#     gli2_equil = []
#     for i,kd in enumerate(KD):
#         print kd
#         x = sim.run(tspan = tspan_equil, initials=initials, param_values= {'kr_tgfbrii_avb3_bind' : kd})
#         gli2_equil.append(x.observables['Gli2_free'][-1])
#         #
#     #     plt.figure()
#     #     plt.plot(tspan_equil, x.observables['Gli2_free'], lw=5, label=r'K$_D$ = %g' % kd)
#     #     plt.legend(loc=0)
#     
#     plt.plot(KD, gli2_equil, lw=5, color=color)
#     
# plt.xscale('log')
# plt.xlabel(r'K$_D$', fontsize=20)
# plt.ylabel(r'[Gli2]$_\mathrm{equil}$', fontsize=20)
# plt.xticks(fontsize=20)
# plt.yticks(fontsize=20)
# plt.tight_layout()
# 
# plt.show()
# quit()
#####

# for i,sp in enumerate(model.species):
#     print i,sp
# quit()

for obs in model.observables:
    plt.figure(obs.name)
    plt.plot(tspan1, x.observables[obs.name], lw=5)  # , label=obs.name)
plt.ylim(ymin=-50, ymax=1250)
plt.xlim(xmax=85)

# save final concentrations
initials = x.species[-1]

# PERTURBATIONS
REMOVE_MECH_FORCES = True
FLUSH_TGFB = False
ADD_1D11 = False
ADD_GANT58 = False

tspan2 = np.append(np.linspace(0, 1, 101), np.linspace(1.1, 40, 130))

# REMOVE MECHANICAL FORCES
if REMOVE_MECH_FORCES:
    x = sim.run(tspan=tspan2, initials=initials,
                param_values={'kf_tgfbrii_avb3_bind': 0})
    #                            'kact_basal_rhoa' : 0.1*kact_basal_rhoa.value})
    #                            'kact_rock' : 0.5*kact_rock.value})
    #  kact_p38mapk
    for obs in model.observables:
        plt.figure(obs.name)
        plt.plot(tspan1[-1]+tspan2, x.observables[obs.name], '-r', lw=5)  # mfc='None', mec='r')

# FLUSH TGFB
if FLUSH_TGFB:
    x = sim.run(tspan=tspan2, initials=initials,
                param_values={'k_tgfb_flush': 1e12})
    for obs in model.observables:
        plt.figure(obs.name)
        plt.plot(tspan1[-1]+tspan2, x.observables[obs.name], '-g', lw=5)  # mfc='None', mec='g')
    
# ADD 1D11 ANTIBODY
if ADD_1D11:
    initials2 = np.array([i for i in initials])
    idx = [str(sp) for sp in model.species].index('x1D11(tgfb=None)')
    initials2[idx] = 100
    ls = ['-', '--', '--', '--']
    for i, Kd in enumerate([0.001, 0.01, 0.1, 1]):
        x = sim.run(tspan=tspan2, initials=initials2,
                    param_values={'kf_1d11_tgfb_bind': 10/Kd,
                                  'kr_1d11_tgfb_bind': 10})
        for obs in model.observables:
            plt.figure(obs.name)
            plt.plot(tspan1[-1]+tspan2, x.observables[obs.name], 'c', ls=ls[i], lw=5)  # mfc='None', mec='c')
    
# ADD GANT58
if ADD_GANT58:
    initials3 = np.array([i for i in initials])
    idx = [str(sp) for sp in model.species].index('GANT58(gli2=None)')
    ls = ['-', '--', '-.', ':']
    for i, g58 in enumerate([1e5, 1e4]):
        initials3[idx] = g58
        x = sim.run(tspan=tspan2, initials=initials3)
        for obs in model.observables:
            plt.figure(obs.name)
            plt.plot(tspan1[-1]+tspan2, x.observables[obs.name], 'm', ls=ls[i], lw=5)  # mfc='None', mec='m')

# Finish plot
for obs in model.observables:
    plt.figure(obs.name)
#     plt.title(obs.name, fontsize=24)
    plt.xlabel('time', fontsize=20)
    plt.ylabel('population', fontsize=20)
    plt.xticks(fontsize=20)
    plt.yticks(fontsize=20)
#     plt.legend(loc=0, fontsize=24)
    plt.tight_layout(pad=3)
    plt.savefig(os.path.join('plots', '%s.pdf' % obs.name), format='pdf')


plt.show()
