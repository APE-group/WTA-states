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

verbose=True
#read general network parameters
crop_pms, plot_pms, sampling_pms = read_crop_and_plot_yaml(is_verbose)

#creating the final set of cropping values by combining with simulation info
crop_pms = set_final_crop_pms(crop_pms, nest_pms, verbose)

default_spectral_analysis_window=plot_pms['spectrogram_window']['default_window']
#switch-case decision of spectral_window_ms based on the value of the crop_duration_ms
if default_spectral_analysis_window:
    spectrogram_window_ms=spectral_window_ms_analysis_switch(crop_pms['duration_ms'])
else:
    spectrogram_window_ms=plot_pms['spectrogram_window']['custom_ms']

#add spectrogram window to sampling_pms
sampling_pms['spectrogram_window_ms']=spectrogram_window_ms

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
if plot_pms['excitatory_rastergram']==True:
    exc_pops_rastegram_plot(cropped_events, num_exc_pop, crop_pms)

#plot of the rastegram of inhibitory neurons
if plot_pms['inhibitory_rastergram']==True:
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

if plot_pms['short_spike_kernel']==True:
    kernel_plot(single_spike_kernel, sampling_rate_Hz)

#collect spike spikes from all populations in a single vector 
single_trace_spike_times = combine_spike_times_in_single_trace(cropped_events, num_exc_pop)

#Calculate combined firing rates
num_exc_pop=nest_pms['network']['num_exc_pop']
num_exc_neu_per_pop=nest_pms['network']['num_exc_neu_per_pop']
time_points, combined_firing_rate = \
    calculate_firing_rate(crop_pms, analysis_pms, single_trace_spike_times, num_exc_pop*num_exc_neu_per_pop)

if plot_pms['instantaneous_firing_rate']==True:
    #plot of combined instantaneous firing rates   
    firing_rates_plot(time_points, combined_firing_rate, crop_pms) 

smoothed_spikes_firing_rate = smooth_signal(combined_firing_rate, single_spike_kernel)
print("len of smoothed signal", len(smoothed_spikes_firing_rate))
print("type of smoothed signal", type(smoothed_spikes_firing_rate))
if plot_pms['smoothed_spikes_firing_rate']==True:
    firing_rates_plot(time_points, smoothed_spikes_firing_rate, crop_pms) 

#dividing the data in segments for reliable spectral analysis
max_plot_freq_Hz=plot_pms['high_freq_spectrum']['max_plot_freq_Hz']
smoothing_length=plot_pms['high_freq_spectrum']['smoothing_lenght'] #frequency samples
if plot_pms['high_freq_spectrum']['plot']==True:
    compute_spectrum_with_error_bands(smoothed_spikes_firing_rate, analysis_pms, max_plot_freq_Hz, smoothing_length)

max_plot_freq_Hz=plot_pms['high_freq_spectrogram']['max_plot_freq_Hz']
#Spectrogram plot
if plot_pms['high_freq_spectrogram']['plot']==True:
    plot_spectrogram(time_points[0],\
                 smoothed_spikes_firing_rate,\
                 analysis_pms, max_plot_freq_Hz)

#smoothing the spike gaussians with synaptic times for ECG/LFP/ECoG like analysis
sampling_rate_Hz = analysis_pms["analysis_sampling_freq_Hz"] #Hz sampling rate
print("sampling_rate_Hz", sampling_rate_Hz)
#duration_ms = 100.0 # Duration of smoothing effect ms
duration_ms = plot_pms['tissue_response_kernel']['duration_ms']
assert(crop_pms["duration_ms"]>=duration_ms)
print("duration_ms", duration_ms)
#std_dev_ms = 12.0  # Standard deviation in ms
std_dev_ms = plot_pms['tissue_response_kernel']['std_dev_ms']
print("std_dev_ms", std_dev_ms)
tissue_response_kernel = gaussian_kernel(duration_ms, sampling_rate_Hz, std_dev_ms)
print("len of kernel",len(tissue_response_kernel))
if plot_pms['tissue_response_kernel']['plot']==True:
    kernel_plot(tissue_response_kernel, sampling_rate_Hz)

#smoothed_spikes_firing_rate = smooth_signal(combined_firing_rate, single_spike_kernel)
tissue_response_rate = smooth_signal(smoothed_spikes_firing_rate, tissue_response_kernel)

print("len of tissue signal", len(tissue_response_rate))
print("type of tissue signal", type(tissue_response_rate))
if plot_pms['tissue_response_rate']==True:
    firing_rates_plot(time_points, tissue_response_rate, crop_pms) 

#dividing the data in segments for reliable spectral analysis
max_plot_freq_Hz=plot_pms['tissue_response_spectrum']['max_plot_freq_Hz']
smoothing_length=plot_pms['tissue_response_spectrum']['smoothing_lenght'] #frequency samples
if plot_pms['tissue_response_spectrum']['plot']==True:
    compute_spectrum_with_error_bands(tissue_response_rate, analysis_pms, max_plot_freq_Hz, smoothing_length)

max_plot_freq_Hz=plot_pms['tissue_response_spectrogram']['max_plot_freq_Hz']
#Spectrogram plot
if plot_pms['tissue_response_spectrogram']['plot']==True:
    plot_spectrogram(time_points[0],\
                 tissue_response_rate,\
                 analysis_pms, max_plot_freq_Hz)

plt.show()
