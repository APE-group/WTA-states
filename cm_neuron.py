#
# Two-compartment neuron creation
#
# First version: 29/04/2024
# Author: Elena Pastorelli, INFN, Rome (IT)
#
# Description: function for the creation of Ca-AdEx
#


from collections import namedtuple

import nest
import numpy as np
import matplotlib.pyplot as plt
import random
import statistics as stat
import sys
import yaml
from dataclasses import dataclass


#sleep = True

@dataclass
class ReceptorIdxs:
    GABA_soma: int = 0
    AMPA_soma: int = 1
    GABA_dist: int = 2
    AMPA_dist: int = 3
    REFRACT_soma: int = 4
    ADAPT_soma: int = 5
    BETA_dist: int = 6
    NMDA_soma: int = 7
    AMPA_NMDA_soma: int = 8
    ALPHAexc_soma: int = 9
    ALPHAinh_soma: int = 10
    ALPHAexc_dist: int = 11
    ALPHAinh_dist: int = 12
    NMDA_dist: int = 13
    AMPA_NMDA_dist: int = 14
    I_soma: int = 15



def init_Ca_AdEx(default_params):

    soma_params = {'C_m': default_params['C_m_s'],                   # [pF] Soma capacitance
                   'g_L': default_params['g_L_s'],                   # [nS] Soma leak conductance
                   'e_L': default_params['e_L_s'],                   # [mV] Soma reversal potential
                   'gbar_Na_Adex': default_params['g_L_s'],          # [nS] Adex conductance
                   'e_Na_Adex': default_params['e_Na_Adex'],         # [mV] Adex threshold
                   'delta_T': default_params['delta_T']              # [mV] Adex slope factor
                   }

    distal_params = {'C_m': default_params['C_m_d'],                 # [pF] Distal capacitance
                     'g_L': default_params['g_L_d'],                 # [nS] Distal leak conductance
                     'g_C': default_params['g_C_d'],                 # [nS] Soma-distal coupling conductance
                     'e_L': default_params['e_L_d'],                 # [mV] Distal reversal potential
                     'gbar_Ca': default_params['gbar_Ca'],           # [nS] Ca maximal conductance
                     'gbar_K_Ca': default_params['gbar_K_Ca'],       # [nS] K_Ca maximal conductance
                     'e_K': default_params['e_K'],                   # [mV] K reversal potential
                     'tau_decay_Ca': default_params['tau_decay_Ca'], # [ms] decay of Ca concentration
                     'phi': default_params['phi'],                   # [-] scale factor
                     'm_half': default_params['m_half'],             # [mV] m half-value for Ca
                     'h_half': default_params['h_half'],             # [mV] h half-value for Ca
                     'm_slope': default_params['m_slope'],           # [-] m slope factor for Ca
                     'h_slope': default_params['h_slope'],           # [-] h slope factor for Ca
                     'tau_m': default_params['tau_m'],               # [ms] m tau decay for Ca
                     'tau_h': default_params['tau_h'],               # [ms] h tau decay dor Ca
                     'tau_m_K_Ca': default_params['tau_m_K_Ca'],     # [ms] m tau decay for K_Ca
                     'Ca_0': default_params['Ca_0'],                 # [mM] Baseline intracellular Ca conc
                     'Ca_th': default_params['Ca_th'],               # [mM] Threshold Ca conc for Ca channel opening
                     'exp_K_Ca': default_params['exp_K_Ca']          # [-] Exponential factor in K_Ca current with Hay dyn
                     }

    other_params = {'d_BAP': default_params['d_BAP'],
                    't_ref': default_params['t_ref'],
                    'V_th': default_params['V_th'],
                    'V_reset': default_params['V_reset'],
                    'tau_w': default_params['tau_w'],
                    'a': default_params['a'],
                    'b': default_params['b'],
                    'e_w': default_params['e_L_s'],
                    'g_w': default_params['g_w'],
                    'w_BAP': default_params['w_BAP']
                    }
        
    return(soma_params, distal_params, other_params)


def init_Ca_AdEx_nestml(default_params, multi_comp=False):
    if multi_comp:
        soma_params = {
            'C_m': default_params['C_m_s'],                   # [pF] Soma capacitance
            'g_L': default_params['g_L_s'],                   # [nS] Soma leak conductance
            'e_L': default_params['e_L_s'],                   # [mV] Soma reversal potential
            'g_spike': default_params['g_L_s'],          # [nS] Adex conductance
            'v_thr': default_params['e_Na_Adex'],         # [mV] Adex threshold
            'delta_thr': default_params['delta_T'],              # [mV] Adex slope factor
            # 'd_BAP': default_params['d_BAP'],             # delay not yet implemented
            'refr_period': default_params['t_ref'],
            'v_reset': default_params['V_reset'],
            'tau_w': default_params['tau_w'],
            'sth_a': default_params['a'],
            'b': default_params['b'],
            'e_adapt': default_params['e_L_s'],
            'g_adapt': default_params['g_w'],
        }

        distal_params = {
            'C_m': default_params['C_m_d'],                 # [pF] Distal capacitance
            'g_L': default_params['g_L_d'],                 # [nS] Distal leak conductance
            'g_C': default_params['g_C_d'],                 # [nS] Soma-distal coupling conductance
            'e_L': default_params['e_L_d'],                 # [mV] Distal reversal potential
            'gbar_Ca': default_params['gbar_Ca'],           # [nS] Ca maximal conductance
            'gbar_K': default_params['gbar_K_Ca'],       # [nS] K_Ca maximal conductance
            'e_K': default_params['e_K'],                   # [mV] K reversal potential
            'tau_Ca': default_params['tau_decay_Ca'], # [ms] decay of Ca concentration
            'phi_Ca': default_params['phi'],                   # [-] scale factor
            'm_half_Ca': default_params['m_half'],             # [mV] m half-value for Ca
            'h_half_Ca': default_params['h_half'],             # [mV] h half-value for Ca
            'm_slope_Ca': default_params['m_slope'],           # [-] m slope factor for Ca
            'h_slope_Ca': default_params['h_slope'],           # [-] h slope factor for Ca
            'tau_m_Ca': default_params['tau_m'],               # [ms] m tau decay for Ca
            'tau_h_Ca': default_params['tau_h'],               # [ms] h tau decay dor Ca
            'tau_m_K': default_params['tau_m_K_Ca'],     # [ms] m tau decay for K_Ca
            'c_Ca0': default_params['Ca_0'],                 # [mM] Baseline intracellular Ca conc
            'c_Ca_thr': default_params['Ca_th'],               # [mM] Threshold Ca conc for Ca channel opening
            'exp_K_Ca': default_params['exp_K_Ca'],
            'w_bap': default_params['w_BAP'],          # [-] Exponential factor in K_Ca current with Hay dyn
        }
    else:
        soma_params = {
            'C_m': default_params['C_m'],                   # [pF] Soma capacitance
            'g_L': default_params['g_L'],                   # [nS] Soma leak conductance
            'e_L': default_params['E_L'],                   # [mV] Soma reversal potential
            'g_spike': default_params['g_L'],          # [nS] Adex conductance
            'v_thr': default_params['V_th'],         # [mV] Adex threshold
            'delta_thr': default_params['Delta_T'],              # [mV] Adex slope factor
            # 'd_BAP': default_params['d_BAP'],             # delay not yet implemented
            'refr_period': default_params['t_ref'],
            'v_reset': default_params['V_reset'],
            'tau_w': default_params['tau_w'],
            'sth_a': default_params['a'],
            'b': default_params['b'],
            'e_adapt': default_params['E_L'],
            'g_adapt': 1., # educated guess
        }
        distal_params = None

    other_params = {'V_th': default_params['V_th']}
        
    return(soma_params, distal_params, other_params)

def create_receptor_mapping(multi_comp=True):        
    if multi_comp:        
        # receptor mapping
        ri = ReceptorIdxs(
            # receptors get assigned an index which corresponds to the order in which they
            # are added. For clearer bookkeeping, we explicitly define these indices here.
            ALPHAexc_soma=0, 
            ALPHAinh_soma=1, 
            ALPHAexc_dist=2, 
            ALPHAinh_dist=3, 
            AMPA_NMDA_dist=4, 
            I_soma=5,
            # other receptors types are not used, assign something that will generate an error when used
            GABA_soma=float('nan'), 
            AMPA_soma=float('nan'), 
            GABA_dist=float('nan'), 
            AMPA_dist=float('nan'), 
            REFRACT_soma=float('nan'), 
            ADAPT_soma=float('nan'), 
            BETA_dist=float('nan'), 
            NMDA_soma=float('nan'), 
            AMPA_NMDA_soma=float('nan'), 
            NMDA_dist=float('nan'), 
        )
    else:
        # receptor mapping
        ri = ReceptorIdxs(
            # receptors get assigned an index which corresponds to the order in which they
            # are added. For clearer bookkeeping, we explicitly define these indices here.
            ALPHAexc_soma=0, 
            ALPHAinh_soma=1, 
            AMPA_NMDA_soma=2, 
            I_soma=3,
            # other receptors types are not used, assign something that will generate an error when used
            GABA_soma=float('nan'), 
            AMPA_soma=float('nan'), 
            GABA_dist=float('nan'), 
            AMPA_dist=float('nan'), 
            REFRACT_soma=float('nan'), 
            ADAPT_soma=float('nan'), 
            BETA_dist=float('nan'), 
            NMDA_soma=float('nan'), 
            ALPHAexc_dist=float('nan'), 
            ALPHAinh_dist=float('nan'), 
            AMPA_NMDA_dist=float('nan'), 
            NMDA_dist=float('nan'), 
        )
    return ri


def create_nestml_neuron(n = 1, params = {}, multi_comp = True):
    soma_params, distal_params, other_params = init_Ca_AdEx_nestml(params, multi_comp=multi_comp)
    # initialize the model with two compartments
    cm = nest.Create("ca_adex_2expsyn_nestml", n)
    cm.V_th = other_params['V_th']

    print(f"? multi_comp: {multi_comp}")

    if multi_comp:
        cm.compartments = [
            {"parent_idx": -1, "params": soma_params},
            {"parent_idx":  0, "params": distal_params}
        ]
        # receptors
        receptors = [
            {"comp_idx": 0, "receptor_type": "syn_2exp", "params": {"tau_r_syn": .173, "tau_d_syn": .227, "e_syn": 0.}}, # ALPHAexc_soma
            {"comp_idx": 0, "receptor_type": "syn_2exp", "params": {"tau_r_syn": 1.73, "tau_d_syn": 2.27, "e_syn": -85.}}, # ALPHAinh_soma
            {"comp_idx": 1, "receptor_type": "syn_2exp", "params": {"tau_r_syn": .173, "tau_d_syn": .227, "e_syn": 0.}}, # ALPHAexc_dist
            {"comp_idx": 1, "receptor_type": "syn_2exp", "params": {"tau_r_syn": 1.73, "tau_d_syn": 2.27, "e_syn": -85.}}, # ALPHAinh_dist
            {"comp_idx": 1, "receptor_type": "ampa_nmda"}, # AMPA_NMDA_dist
            {"comp_idx": 0, "receptor_type": "inp"}, # I_SOMA
        ]
        cm.receptors = receptors

    else:
        cm.compartments =[
            {"parent_idx": -1, "params": soma_params},
        ]
        # receptors
        receptors = [
            {"comp_idx": 0, "receptor_type": "syn_2exp", "params": {"tau_r_syn": .173, "tau_d_syn": .227, "e_syn": 0.}}, # ALPHAexc_soma
            {"comp_idx": 0, "receptor_type": "syn_2exp", "params": {"tau_r_syn": 1.73, "tau_d_syn": 2.27, "e_syn": -85.}}, # ALPHAinh_soma
            {"comp_idx": 0, "receptor_type": "ampa_nmda"}, # AMPA_NMDA_soma
            {"comp_idx": 0, "receptor_type": "inp"}, # I_SOMA
        ]
        cm.receptors = receptors

    return cm


def create_cm_neuron(n = 1, params = {}):
    soma_params, distal_params, other_params = init_Ca_AdEx(params)
    print('soma params: ', soma_params)
    print('distal params: ', distal_params)
    print('other params: ', other_params)
    # create a model with two compartments
    print("Creating %d Ca-AdEx neurons..." %n)
    cm = nest.Create('cm_default', n)
    cm.compartments = [
        {"parent_idx": -1, "params": soma_params},
        {"parent_idx":  0, "params": distal_params}
    ]

    V_th = other_params['V_th']
    V_reset = other_params['V_reset']
    t_ref = other_params['t_ref']
    tau_w = other_params['tau_w']
    a = other_params['a']
    b = other_params['b']
    e_w = other_params['e_w']
    g_w = other_params['g_w']

    # spike threshold
    nest.SetStatus(cm, {'V_th': V_th})   # Spike threshold (mV)
    nest.SetStatus(cm, {'V_reset': V_reset})   # Voltage reset (mV)

    # receptors
    receptors = [
        {"comp_idx": 0, "receptor_type": "GABA", "params": {"tau_r_GABA": .2, "tau_d_GABA": 2., "e_GABA": -85.}}, # GABA_soma,
        {"comp_idx": 0, "receptor_type": "AMPA", "params": {"tau_r_AMPA": .2, "tau_d_AMPA": 3., "e_AMPA": 0.}}, # AMPA_soma
        {"comp_idx": 1, "receptor_type": "GABA", "params": {"tau_r_GABA": .2, "tau_d_GABA": 2., "e_GABA": -85.}}, # GABA_dist
        {"comp_idx": 1, "receptor_type": "AMPA", "params": {"tau_r_AMPA": .2, "tau_d_AMPA": 3., "e_AMPA": 0.}}, # AMPA_dist
        {"comp_idx": 0, "receptor_type": "REFRACT", "params": {"t_ref": t_ref, "g_refract": 100000., "V_reset": V_reset}}, # REFRACT_soma
        {"comp_idx": 0, "receptor_type": "ADAPT", "params": {"tau_w": tau_w, "a": a, "b": b, "e_w": e_w, "g_w": g_w}}, #  ADAPT_soma
        {"comp_idx": 1, "receptor_type": "BETA", "params": {"tau_r_BETA": 2., "tau_d_BETA": 5.}}, # BETA_dist
        {"comp_idx": 0, "receptor_type": "NMDA"}, # NMDA_soma
        {"comp_idx": 0, "receptor_type": "AMPA_NMDA"}, # AMPA_NMDA_soma
        {"comp_idx": 0, "receptor_type": "AMPA", "params": {"tau_r_AMPA": .173, "tau_d_AMPA": .227, "e_AMPA": 0.}}, # ALPHAexc_soma
        {"comp_idx": 0, "receptor_type": "AMPA", "params": {"tau_r_AMPA": 1.73, "tau_d_AMPA": 2.27, "e_AMPA": -85.}}, # ALPHAinh_soma
        {"comp_idx": 1, "receptor_type": "AMPA", "params": {"tau_r_AMPA": .173, "tau_d_AMPA": .227, "e_AMPA": 0.}}, # ALPHAexc_dist
        {"comp_idx": 1, "receptor_type": "AMPA", "params": {"tau_r_AMPA": 1.73, "tau_d_AMPA": 2.27, "e_AMPA": -85.}}, # ALPHAinh_dist
        {"comp_idx": 1, "receptor_type": "NMDA"}, #  NMDA_dist
        {"comp_idx": 1, "receptor_type": "AMPA_NMDA"}, # AMPA_NMDA_dist
    ]
    cm.receptors = receptors

    # receptors get assigned an index which corresponds to the order in which they
    # are added. For clearer bookkeeping, we explicitly define these indices here.
    GABA_soma, AMPA_soma, GABA_dist, AMPA_dist, REFRACT_soma, ADAPT_soma, BETA_dist, NMDA_soma, AMPA_NMDA_soma, ALPHAexc_soma, ALPHAinh_soma, ALPHAexc_dist, ALPHAinh_dist, NMDA_dist, AMPA_NMDA_dist = 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14

    # connect soma to soma to generate refractoriness
    syn_dict_REFRACT = {
        'synapse_model': 'REFRACT_syn',
        'weight': 1., 
        'delay': .1, 
        'receptor_type': REFRACT_soma}
    for i in range(n):
        nest.Connect(cm[i], cm[i], syn_spec=syn_dict_REFRACT)

    # connect soma to soma to reproduce adaptation currents
    syn_dict_ADAPT = {
        'synapse_model': 'ADAPT_syn',
        'weight': 1., 
        'delay': .1, 
        'receptor_type': ADAPT_soma}
    for i in range(n):
        nest.Connect(cm[i], cm[i], syn_spec=syn_dict_ADAPT)

    # connect soma to distal with delay to reproduce back propagation currents
    w_BAP = other_params['w_BAP']
    d_BAP = other_params['d_BAP']
    syn_dict_BAP = {
        'synapse_model': 'BAP_syn',
        'weight': w_BAP, 
        'delay': d_BAP, 
        'receptor_type': AMPA_dist}
    for i in range(n):
        nest.Connect(cm[i], cm[i], syn_spec=syn_dict_BAP)

    return(cm)
