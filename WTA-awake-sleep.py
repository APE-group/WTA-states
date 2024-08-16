#!/usr/bin/env python
# coding: utf-8

#imports
import nest
import yaml
import nest.raster_plot
import numpy as np
import matplotlib.pyplot as plt
from scipy.signal import spectrogram, butter, filtfilt
from analysis_functions import *
from asserting_functions import *
from yaml_io import *
from prepare_nest_parameters import *
#from spectrum_with_bands import *
from nest_reset_create_connect_simulate import *

is_verbose = False
debug_mode = True

#total sim, resolution and recording times
times = read_sim_and_recording_times_yaml(is_verbose)

#read general network parameters
config = read_general_config_yaml(is_verbose)

#prepare all simulation parameters
nest_pms={}
nest_pms = nest_parameters_preparation(times, config, is_verbose, nest_pms)
print("nest_pms",nest_pms)

NEST_version = nest.__version__
if NEST_version == "3.7.0" and nest_pms["use_single_compartment_environment"]==False:
    print("ASSERTION ERROR: Ca-AdEx multi-compartment neuron not supported by this NEST version", NEST_version)
    assert(False)

num_threads=4
sim_completed, spike_recorders, inh_spike_recorder = nest_reset_create_connect_simulate(nest_pms,num_threads)
print("sim_completed", sim_completed)

d_inh = nest.GetStatus(inh_spike_recorder, "events")[0]


#before analysis, preliminary sim look 
preliminary_sim_look(debug_mode,nest_pms, spike_recorders, inh_spike_recorder, nest_pms["recording_pms"])   

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

#plot of the rastegram of excitatory neurons
num_exc_pop=nest_pms['network']['num_exc_pop']
exc_pops_rastegram_plot(cropped_events, num_exc_pop, crop_pms)

#plot of the rastegram of inhibitory neurons
inh_rastergram_plot(cropped_inh_events, crop_pms)

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

#collect spike spikes from all populations in a single vector 
single_trace_spike_times = combine_spike_times_in_single_trace(cropped_events, num_exc_pop)

#Calculate combined firing rates
num_exc_pop=nest_pms['network']['num_exc_pop']
num_exc_neu_per_pop=nest_pms['network']['num_exc_neu_per_pop']
time_points, combined_firing_rate = \
    calculate_firing_rate(crop_pms, analysis_pms, single_trace_spike_times, num_exc_pop*num_exc_neu_per_pop)


instantaneous_firing_rate_plot=False
if instantaneous_firing_rate_plot:
    #plot of combined instantaneous firing rates   
    firing_rates_plot(time_points, combined_firing_rate, crop_pms) 

smoothed_spikes_firing_rate = smooth_signal(combined_firing_rate, single_spike_kernel)
print("len of smoothed signal", len(smoothed_spikes_firing_rate))
print("type of smoothed signal", type(smoothed_spikes_firing_rate))
firing_rates_plot(time_points, smoothed_spikes_firing_rate, crop_pms) 

#dividing the data in segments for reliable spectral analysis
max_plot_freq_Hz=250.0
smoothing_length=20 #frequency samples
compute_spectrum_with_error_bands(smoothed_spikes_firing_rate, analysis_pms, max_plot_freq_Hz, smoothing_length)

max_plot_freq_Hz=100

#Spectrogram plot
plot_spectrogram(time_points[0],\
                 smoothed_spikes_firing_rate,\
                 analysis_pms, max_plot_freq_Hz)

#Spectrogram plot
#plot_spectrogram(time_points[0], combined_firing_rate, analysis_pms)

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

#smoothed_spikes_firing_rate = smooth_signal(combined_firing_rate, single_spike_kernel)
tissue_response_rate = smooth_signal(smoothed_spikes_firing_rate, tissue_response_kernel)

print("len of tissue signal", len(tissue_response_rate))
print("type of tissue signal", type(tissue_response_rate))
firing_rates_plot(time_points, tissue_response_rate, crop_pms) 


#dividing the data in segments for reliable spectral analysis
max_plot_freq_Hz=40
smoothing_length=4 #frequency samples
compute_spectrum_with_error_bands(tissue_response_rate, analysis_pms, max_plot_freq_Hz, smoothing_length)

max_plot_freq_Hz=40
#Spectrogram plot
plot_spectrogram(time_points[0],\
                 tissue_response_rate,\
                 analysis_pms, max_plot_freq_Hz)

#Spectrogram plot
#plot_spectrogram(time_points[0], combined_firing_rate, analysis_pms)

plt.show()






