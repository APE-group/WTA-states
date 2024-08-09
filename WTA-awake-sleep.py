#!/usr/bin/env python
# coding: utf-8

# In[1]:


#imports
import nest
import yaml
import nest.raster_plot
import numpy as np
import matplotlib.pyplot as plt
from scipy.signal import spectrogram, butter, filtfilt
from analysis_functions import *
from asserting_functions import *
#from spectrum_with_bands import *
from nest_reset_create_connect_simulate import *
from yaml_io import *


# In[2]:


NEST_version = nest.__version__
is_verbose=True
debug_mode = True
high_noise = True
nest_pms={}

#total sim, resolution and recording times
times = read_sim_and_recording_times_yaml(is_verbose)

#read general network parameters
config = read_general_config_yaml(is_verbose)

# Extract config subsets of params 
brain_state = config['brain_state']
use_single_compartment_environment = config['use_single_compartment_environment']
exc_pms_file_name = config['exc_neu_params_filename']
inh_pms_file_name = config['inh_neu_params_filename']
network = config['network']
cf = config['cf']
weights = config['weights']
poisson = config['poisson']
dc_exc = config['dc_exc']
dc_inh = config['dc_inh']

print("NEST version",NEST_version)


# In[3]:


#reading neural parameters

exc_pms = read_neural_parameters(exc_pms_file_name, is_verbose)
inh_pms = read_neural_parameters(inh_pms_file_name, is_verbose)

nest_pms["exc_pms"]=inh_pms
nest_pms["exc_pms"]=exc_pms
exc_t_ref_ms=0.0


# In[4]:


num_exc_pop=network['num_exc_pop']
num_exc_neu_per_pop=network['num_exc_neu_per_pop']
# Calculate number of inhibitory neurons based on excitatory populations and neurons per pop
network['num_inh_neu'] = int(network['num_exc_pop'] * network['num_exc_neu_per_pop'] / 4)
num_inh_neu = network['num_inh_neu']

#setting a minimum sinaptic delay if exc_t_ref_ms<=2.0
min_syn_delay_ms = 3.0 if exc_t_ref_ms else 0.0


# In[5]:


use_recurrency = network["use_exc_recurrency"]
if is_verbose:
    print("use_exc_recurrency:",use_recurrency)
    print("weights:",weights)


# In[6]:


if use_recurrency:
    if brain_state == "awake":
        weights=weights['awake']
    else:
        weights=weights['NREM'] 
    recurrent_weight = weights['recurrent_weight'] * cf if use_recurrency else 0.0
    inh_to_exc_weight = weights["inh_to_exc_weight"] * cf
    exc_to_inh_weight = weights["exc_to_inh_weight"]
    inh_to_inh_weight = weights["inh_to_inh_weight"]


# In[7]:


if is_verbose:
    print("recurrent_weight:", recurrent_weight)
    print("inh_to_exc_weight", inh_to_exc_weight)
    print("exc_to_inh_weight", exc_to_inh_weight)
    print("inh_to_inh_weight", inh_to_inh_weight)


# In[8]:


use_poisson_generators = network["use_poisson_generators"]
if use_poisson_generators:
    num_poisson_generators=network["num_poisson_generators"]
    # Set stimulation parameters based on brain state
    if brain_state == "awake":
        poisson=poisson['awake']
    else:
        if brain_state == "NREM":
            poisson=poisson['NREM']
        else:
            assert(false==true)
    poisson_spreading_factor=poisson['poisson_noise_spreading_factor']
    poisson_rate = poisson['poisson_basic_rate']/ poisson_spreading_factor
    poisson_weight = poisson['poisson_weight'] * poisson_spreading_factor * cf
    poisson_delta = poisson["poisson_delta"]

# Create rate arrays for Poisson generators
#poisson_rate_arrays = [poisson_rate + i * poisson_delta for i in range(network["num_exc_neu_per_pop"])]
#print("number of poisson generators",len(poisson_rate_arrays))


# In[9]:


#dc inhjectors to excitatory neurons
use_dc_exc_injectors = network["use_dc_exc_injectors"]
if use_dc_exc_injectors: 
    dc_exc_amplitudes = dc["amplitudes"]
    dc_exc_weight = dc["weight"]
    dc_exc_delay = dc["delay"]
    dc_exc_start = 0.0
    dc_exc_stop = times["sim_pms"]["stop_ms"]

#dc inhjectors to inhibitory neurons
use_dc_inh_injector = network["use_dc_inh_injector"]  # DC injection for inhibitory neurons
if use_dc_inh_injector:
    dc_inh_amplitude = dc_inh["amplitude"] # DC amplitude for inhibitory neurons
    dc_inh_weight = dc_inh["weight"]
    dc_inh_delay = dc_inh["delay"]
    dc_inh_start = 0.0
    dc_inh_stop = times["sim_pms"]["stop_ms"]


# In[10]:


if is_verbose:
    print ("exc neu params before brain-state specific tuning:", exc_pms)


# In[11]:


if is_verbose:
    print ("inh neu params before brain-state specific tuning:", inh_pms)


# In[12]:


assert(exc_pms['neuron_kind']=="excitatory")
if exc_pms['multi_compartment']==False:
    print("setting brain state INDEPENDENT neural params in single-comp neu")
    # Set brain state independent parameters for single-compartment neuron
    exc_neu_params={}
    exc_neu_params['receptors']={} 
    exc_neu_params['equation_params'] = {
        "a": exc_pms['equation_params']['a'],
        "t_ref": exc_pms['equation_params']['t_ref'],
        "Delta_T": exc_pms['equation_params']['Delta_T'],
        "C_m": exc_pms['equation_params']['C_m'],
        "g_L": exc_pms['equation_params']['g_L'],
        "tau_w": exc_pms['equation_params']['tau_w'],
        "V_th": exc_pms['equation_params']['V_th'],
        "V_peak": exc_pms['equation_params']['V_peak'],
    }
    print("setting brain state DEPENDENT neural params in single-comp neu")
    # Check if brain_state variable is set correctly
    assert(brain_state in ["awake", "NREM"])
    # Set brain state dependent parameters for single-compartment neuron
    if brain_state == "awake":
        brain_state_dependent_params = {
            "b": exc_pms['equation_params'].get('b_awake', None),
            "E_L": exc_pms['equation_params'].get('E_L_awake', None),
            "V_reset": exc_pms['equation_params'].get('V_reset_awake', None)
        }
        print(f"Brain state 'awake' parameters: {brain_state_dependent_params}")
        exc_neu_params['equation_params'].update(brain_state_dependent_params)
    elif brain_state == "NREM":
        brain_state_dependent_params = {
            "b": exc_pms['equation_params'].get('b_NREM', None),
            "E_L": exc_pms['equation_params'].get('E_L_NREM', None),
            "V_reset": exc_pms['equation_params'].get('V_reset_NREM', None)
        }
        exc_neu_params['equation_params'].update(brain_state_dependent_params)
else:
    assert(exc_pms['multi_compartment']==True)
    assert(use_single_compartment_environment==False)
    if is_verbose:
        print("setting brain state dependent neural params in multi-comp neu")
    exc_neu_params = {}
    exc_neu_params['receptors']=exc_pms['receptors']
    exc_neu_params['equation_params']={}
    exc_neu_params['equation_params'].update(exc_pms['equation_params'])
    # Set brain state dependent parameters for multi-compartment neuron
    if brain_state == "awake":
        print("for multi-compartment neurons the default set of parameters is the awake one")
    elif brain_state == "NREM":
        brain_state_dependent_params = {
            'b': exc_pms['NREM_changers'].get('b_NREM', None),
            'e_L_s': exc_pms['equation_params'].get('e_L_s', None) + exc_pms['NREM_changers'].get('delta_e_L_s_NREM', 0),
            'e_L_d': exc_pms['equation_params'].get('e_L_d', None) + exc_pms['NREM_changers'].get('delta_e_L_d_NREM', 0),
            'g_C_d': exc_pms['equation_params'].get('g_C_d', None) * exc_pms['NREM_changers'].get('g_C_d_NREM_multiplier', 1)
        }
        print(f"Brain state 'sleep' parameters for multi-compartment neuron: {brain_state_dependent_params}")
        exc_neu_params['equation_params'].update(brain_state_dependent_params)

if is_verbose:
    print("exc neu params AFTER brain-state specific tuning:", exc_neu_params)


# In[13]:


assert(inh_pms['neuron_kind']=="inhibitory")
assert(inh_pms['multi_compartment']==False)
inh_neu_params=inh_pms['equation_params']


# In[14]:


if is_verbose:
    print ("inh neu params AFTER brain-state specific tuning:", inh_neu_params)


# In[15]:


nest_pms={}
# Add contextual poisson signal configuration to nest_pms
if 'contextual_poisson' in config:
    nest_pms['contextual_poisson'] = config['contextual_poisson']
nest_pms["sim_pms"]=times["sim_pms"]
nest_pms["use_single_compartment_environment"]= use_single_compartment_environment
nest_pms["exc_neu_params"]=exc_neu_params
nest_pms["inh_neu_params"]=inh_neu_params
nest_pms["exc_t_ref_ms"]=exc_t_ref_ms
nest_pms["recording_pms"]=times["recording_pms"]

#how many neurons
nest_pms["network"]={}
nest_pms["network"]["num_exc_neu_per_pop"]=num_exc_neu_per_pop
nest_pms["network"]["num_exc_pop"]=num_exc_pop
nest_pms["network"]["num_inh_neu"]=num_inh_neu

#synaptic weights
nest_pms["network"]["use_exc_recurrency"]=use_recurrency
nest_pms["network"]["weights"]={}
nest_pms["network"]["weights"]["exc_to_exc"] = recurrent_weight 
nest_pms["network"]["weights"]["inh_to_exc_weight"] = inh_to_exc_weight
nest_pms["network"]["weights"]["exc_to_inh_weight"] = exc_to_inh_weight
nest_pms["network"]["weights"]["inh_to_inh_weight"] = inh_to_inh_weight
nest_pms["network"]["min_syn_delay_ms"]=min_syn_delay_ms
#poisson generators
nest_pms["network"]["use_poisson_generators"]=use_poisson_generators
if(use_poisson_generators):
    nest_pms["network"]["num_poisson_generators"]=num_poisson_generators
    nest_pms["poisson"]={}
    nest_pms["poisson"]["poisson_noise_spreading_factor"]=poisson_spreading_factor
    nest_pms["poisson"]["poisson_rate"]=poisson_rate
    nest_pms["poisson"]["poisson_weight"]=poisson_weight
    nest_pms["poisson"]["poisson_delta"]=poisson_delta

#DC exc injectors
nest_pms["use_dc_exc_injectors"] = use_dc_exc_injectors
if(use_dc_exc_injectors):
    nest_pms["dc_exc"]={}
    nest_pms["dc_exc"]["start_ms"] = dc_exc_start
    nest_pms["dc_exc"]["stop_ms"] = dc_exc_stop
    nest_pms["dc_exc"]["dc_exc_amplitudes"] = dc_exc_amplitudes
    nest_pms["dc_exc"]["dc_exc_weight"] = dc_exc_weight
    nest_pms["dc_exc"]["dc_exc_delay_ms"] = dc_exc_delay

#DC inh injectors
nest_pms["use_dc_inh_injector"] = use_dc_inh_injector
if use_dc_inh_injector:
    nest_pms["dc_inh"]={}
    nest_pms["dc_inh"]["start_ms"]=dc_inh_start 
    nest_pms["dc_inh"]["stop_ms"]=dc_inh_stop 
    nest_pms["dc_inh"]["dc_inh_amplitude"] = dc_inh_amplitude
    nest_pms["dc_inh"]["weight"]=dc_inh_weight
    nest_pms["dc_inh"]["delay_ms"]=dc_inh_delay 


# In[16]:


print("nest_pms",nest_pms)


# In[17]:


NEST_version = nest.__version__
if NEST_version == "3.7.0" and nest_pms["use_single_compartment_environment"]==False:
    print("ASSERTION ERROR: Ca-AdEx multi-compartment neuron not supported by this NEST version", NEST_version)
    assert(False)


# In[18]:


num_threads=4
sim_completed, spike_recorders, inh_spike_recorder = nest_reset_create_connect_simulate(nest_pms,num_threads)
print("sim_completed", sim_completed)


# In[19]:


d_inh = nest.GetStatus(inh_spike_recorder, "events")[0]


# In[ ]:





# In[20]:


#before analysis, preliminary sim look 
preliminary_sim_look(debug_mode,nest_pms, spike_recorders, inh_spike_recorder, nest_pms["recording_pms"])   


# In[ ]:





# In[21]:


#setting the analysis parameters
default_cropping=True

if default_cropping:
    crop_start_ms=nest_pms["recording_pms"]["start_ms"] 
    crop_stop_ms=nest_pms["recording_pms"]["stop_ms"] 
else:
    crop_start_ms = 000.0 + nest_pms["recording_pms"]["start_ms"]
    crop_stop_ms= 200.0 + nest_pms["recording_pms"]["start_ms"]

crop_pms={"start_ms":crop_start_ms, 
          "stop_ms":crop_stop_ms}
crop_pms["duration_ms"]=crop_pms["stop_ms"]-crop_pms["start_ms"]

#switch-case decision of spectral_window_ms based on the value of the crop_duration_ms
def spectral_window_ms_analysis_switch(crop_duration_ms):
    if crop_duration_ms<4000.0:
        spectral_window_ms=crop_duration_ms/4.0
        return spectral_window_ms
    if crop_duration_ms <= 8000.0:
        spectral_window_ms = 1000.0
        return spectral_window_ms
    if crop_duration_ms <= 16000.0:
        spectral_window_ms = 2000.0
        return spectral_window_ms
    else:
        spectral_window_ms = 4000.0
        return spectral_window_ms
        
default_spectral_analysis=True

if default_spectral_analysis:
    spectrogram_window_ms=spectral_window_ms_analysis_switch(crop_pms["duration_ms"])
else:
    spectrogram_window_ms=1000.0
   
sampling_pms={"spikes_sampling_window_ms": 0.2,
              "spectrogram_window_ms":spectrogram_window_ms,
              "low_freq_sampling_window_ms":12.5}

analysis_pms={}                                 
analysis_pms.update(derive_analysis_pms_from_sampling_pms(sampling_pms))
  
#assertions on crop, recording and spectrogram params
check_analysis_pms(crop_pms,nest_pms["recording_pms"],analysis_pms)

# Assuming spike_recorders is a list of spike recorder IDs previously created in your NEST simulation
cropped_events = crop_events_from_spike_recorders(crop_pms, spike_recorders)
"""
print("type of cropped_events",type(cropped_events), len(cropped_events))
print(cropped_events)
"""
cropped_inh_events = crop_inh_events(crop_pms, inh_spike_recorder)
"""
print("type of cropped_inh_events",type(cropped_inh_events), len(cropped_inh_events))
print(cropped_inh_events)
"""
#print the set of analysis parameters
analysis_pms_print(nest_pms["recording_pms"],crop_pms,analysis_pms)


# In[22]:


#plot of the rastegram of excitatory neurons
exc_pops_rastegram_plot(cropped_events, num_exc_pop, crop_pms)


# In[23]:


#plot of the rastegram of inhibitory neurons
inh_rastergram_plot(cropped_inh_events, crop_pms)


# In[24]:


#preparing to smooth the single spike with the form of a spike event
sampling_rate_Hz = analysis_pms["analysis_sampling_freq_Hz"] #Hz sampling rate
print("sampling_rate_Hz", sampling_rate_Hz)
duration_ms = 4.0 # Duration of smoothing effect ms
print("duration_ms", duration_ms)
std_dev_ms = 0.5  # Standard deviation in ms
print("std_dev_ms", std_dev_ms)

single_spike_kernel = gaussian_kernel(duration_ms, sampling_rate_Hz, std_dev_ms)
print("len of kernel",len(single_spike_kernel))

kernel_plot(single_spike_kernel, sampling_rate_Hz)


# In[25]:


#collect spike spikes from all populations in a single vector 
single_trace_spike_times = combine_spike_times_in_single_trace(cropped_events, num_exc_pop)


# In[26]:


#Calculate combined firing rates
time_points, combined_firing_rate = \
    calculate_firing_rate(crop_pms, analysis_pms, single_trace_spike_times, num_exc_pop*num_exc_neu_per_pop)


# In[27]:


instantaneous_firing_rate_plot=False
if instantaneous_firing_rate_plot:
    #plot of combined instantaneous firing rates   
    firing_rates_plot(time_points, combined_firing_rate, crop_pms) 


# In[28]:


smoothed_spikes_firing_rate = smooth_signal(combined_firing_rate, single_spike_kernel)
print("len of smoothed signal", len(smoothed_spikes_firing_rate))
print("type of smoothed signal", type(smoothed_spikes_firing_rate))
firing_rates_plot(time_points, smoothed_spikes_firing_rate, crop_pms) 


# In[29]:


#dividing the data in segments for reliable spectral analysis
max_plot_freq_Hz=250.0
smoothing_length=20 #frequency samples
compute_spectrum_with_error_bands(smoothed_spikes_firing_rate, analysis_pms, max_plot_freq_Hz, smoothing_length)


# In[30]:


max_plot_freq_Hz=100

#Spectrogram plot
plot_spectrogram(time_points[0],\
                 smoothed_spikes_firing_rate,\
                 analysis_pms, max_plot_freq_Hz)

#Spectrogram plot
#plot_spectrogram(time_points[0], combined_firing_rate, analysis_pms)


# In[31]:


#smoothing the spike gaussians with synaptic times for ECG/LFP/ECoG like analysis
sampling_rate_Hz = analysis_pms["analysis_sampling_freq_Hz"] #Hz sampling rate
print("sampling_rate_Hz", sampling_rate_Hz)
duration_ms = 100.0 # Duration of smoothing effect ms
assert(crop_pms["duration_ms"]>=duration_ms)
print("duration_ms", duration_ms)
std_dev_ms = 12.0  # Standard deviation in ms
print("std_dev_ms", std_dev_ms)
tissue_response_kernel = gaussian_kernel(duration_ms, sampling_rate_Hz, std_dev_ms)
print("len of kernel",len(tissue_response_kernel))

kernel_plot(tissue_response_kernel, sampling_rate_Hz)


# In[32]:


#smoothed_spikes_firing_rate = smooth_signal(combined_firing_rate, single_spike_kernel)
tissue_response_rate = smooth_signal(smoothed_spikes_firing_rate, tissue_response_kernel)


# In[33]:


print("len of tissue signal", len(tissue_response_rate))
print("type of tissue signal", type(tissue_response_rate))
firing_rates_plot(time_points, tissue_response_rate, crop_pms) 


# In[34]:


#dividing the data in segments for reliable spectral analysis
max_plot_freq_Hz=40
smoothing_length=4 #frequency samples
compute_spectrum_with_error_bands(tissue_response_rate, analysis_pms, max_plot_freq_Hz, smoothing_length)


# In[35]:


max_plot_freq_Hz=40
#Spectrogram plot
plot_spectrogram(time_points[0],\
                 tissue_response_rate,\
                 analysis_pms, max_plot_freq_Hz)

#Spectrogram plot
#plot_spectrogram(time_points[0], combined_firing_rate, analysis_pms)


# In[36]:


plt.show()


# In[ ]:





# In[ ]:





# In[ ]:




