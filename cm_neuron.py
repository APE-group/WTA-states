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


#sleep = True


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
                    'w_BAP': default_params['w_BAP']
                    }
        
    return(soma_params, distal_params, other_params)

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

    # spike threshold
    nest.SetStatus(cm, {'V_th': V_th})   # Spike threshold (mV)
    nest.SetStatus(cm, {'V_reset': V_reset})   # Voltage reset (mV)

    # receptors
    receptors = [
        {"comp_idx": 0, "receptor_type": "GABA", "params": {"tau_r_GABA": .2, "tau_d_GABA": 2., "e_GABA": -85.}},
        {"comp_idx": 0, "receptor_type": "AMPA", "params": {"tau_r_AMPA": .2, "tau_d_AMPA": 3., "e_AMPA": 0.}},
        {"comp_idx": 1, "receptor_type": "GABA", "params": {"tau_r_GABA": .2, "tau_d_GABA": 2., "e_GABA": -85.}},
        {"comp_idx": 1, "receptor_type": "AMPA", "params": {"tau_r_AMPA": .2, "tau_d_AMPA": 3., "e_AMPA": 0.}},
        {"comp_idx": 0, "receptor_type": "REFRACT", "params": {"t_ref": t_ref, "g_refract": 100000., "V_reset": V_reset}},
        {"comp_idx": 0, "receptor_type": "ADAPT", "params": {"tau_w": tau_w, "a": a, "b": b, "e_w": e_w}},
        {"comp_idx": 1, "receptor_type": "BETA", "params": {"tau_r_BETA": 2., "tau_d_BETA": 5.}},
        {"comp_idx": 0, "receptor_type": "NMDA"},
        {"comp_idx": 0, "receptor_type": "AMPA_NMDA"},
        {"comp_idx": 0, "receptor_type": "AMPA", "params": {"tau_r_AMPA": .173, "tau_d_AMPA": .227, "e_AMPA": 0.}},
        {"comp_idx": 0, "receptor_type": "AMPA", "params": {"tau_r_AMPA": 1.73, "tau_d_AMPA": 2.27, "e_AMPA": -85.}},
    ]
    
    cm.receptors = receptors

    # receptors get assigned an index which corresponds to the order in which they
    # are added. For clearer bookkeeping, we explicitly define these indices here.
    GABA_soma, AMPA_soma, GABA_dist, AMPA_dist, REFRACT_soma, ADAPT_soma, BETA_dist, NMDA, AMPA_NMDA, ALPHAexc_soma, ALPHAinh_soma = 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10

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
