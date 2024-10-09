#!/usr/bin/env python
# coding: utf-8
import numpy as np
import matplotlib.pyplot as plt
from yaml_io import *
from asserting_functions import *
from analysis_functions import *

def prepare_crop_plot_sampling_activityAnalysis_parameters(directories_and_list_of_yamls, nest_pms, verbose):
    # Read general network parameters
    crop_pms, plot_pms, sampling_pms = read_crop_and_plot_yaml(directories_and_list_of_yamls, verbose)
    
    # Create the final set of cropping values by combining with simulation info
    crop_pms = set_final_crop_pms(crop_pms, nest_pms, verbose)
    
    default_spectral_analysis_window = plot_pms['spectrogram_window']['default_window']
    # Decide spectral_window_ms based on the value of the crop_duration_ms
    if default_spectral_analysis_window:
        spectrogram_window_ms = spectral_window_ms_analysis_switch(crop_pms['duration_ms'])
    else:
        spectrogram_window_ms = plot_pms['spectrogram_window']['custom_ms']
    
    # Add spectrogram window to sampling_pms
    sampling_pms['spectrogram_window_ms'] = spectrogram_window_ms
    
    analysis_pms = {}                                 
    analysis_pms.update(derive_analysis_pms_from_sampling_pms(sampling_pms))
      
    # Assertions on crop, recording, and spectrogram params
    check_analysis_pms(crop_pms, nest_pms["recording_pms"], analysis_pms)
    
    # Print the set of analysis parameters
    analysis_pms_print(nest_pms["recording_pms"], crop_pms, analysis_pms)

    return crop_pms, plot_pms, sampling_pms, analysis_pms

def produce_rastegrams(nest_pms, plot_pms, crop_pms, cropped_events, cropped_inh_events, verbose):
    # Plot of the rastergram of excitatory neurons
    num_exc_pop = nest_pms['network']['num_exc_pop']
    if plot_pms['excitatory_rastergram']:
        exc_pops_rastegram_plot(cropped_events, num_exc_pop, crop_pms, plot_pms)
    # Plot of the rastergram of inhibitory neurons
    if plot_pms['inhibitory_rastergram']:
        inh_rastergram_plot(cropped_inh_events, crop_pms, plot_pms)
    return

def prepare_single_spike_kernel(analysis_pms, plot_pms, crop_pms, verbose):
    # Creating the kernel to smooth an instantaneous spike time into a spike-like waveform
    sampling_rate_Hz = analysis_pms["analysis_sampling_freq_Hz"]  # Hz sampling rate
    if verbose:
        print("in prepare_single_spike_kernel: sampling_rate_Hz", sampling_rate_Hz)
    
    duration_ms = plot_pms['spike_activation_kernel']['duration_ms']
    if crop_pms["duration_ms"] < duration_ms:
        print(f"Warning: Cropping duration ({crop_pms['duration_ms']} ms) is shorter than the spike activation kernel duration ({duration_ms} ms). Cannot perform spike-like waveform analysis.")
        return None
    
    if verbose:
        print("in prepare_single_spike_kernel: duration_ms", duration_ms)
        
    std_dev_ms = plot_pms['spike_activation_kernel']['std_dev_ms']  # Standard deviation in ms
    if verbose:
        print("in prepare_single_spike_kernel: std_dev_ms", std_dev_ms)
    # Create the kernel
    single_spike_kernel = gaussian_kernel(duration_ms, sampling_rate_Hz, std_dev_ms)
    if verbose:
        print("in prepare_single_spike_kernel: len of kernel", len(single_spike_kernel))
    
    if plot_pms['spike_activation_kernel']['plot']:
        kernel_plot(single_spike_kernel, sampling_rate_Hz, plot_pms, 'spike-like_waveform_')
    return single_spike_kernel

def prepare_tissue_response_kernel(analysis_pms, plot_pms, crop_pms, verbose):
    # Creating a kernel to smooth the spike-like waveforms with tissue- and synaptic-like times
    sampling_rate_Hz = analysis_pms["analysis_sampling_freq_Hz"]  # Hz sampling rate
    if verbose:
        print("in prepare_tissue_response_kernel: sampling_rate_Hz", sampling_rate_Hz)
    
    duration_ms = plot_pms['tissue_response_kernel']['duration_ms']
    if crop_pms["duration_ms"] < duration_ms:
        print(f"Warning: Cropping duration ({crop_pms['duration_ms']} ms) is shorter than the tissue response kernel duration ({duration_ms} ms). Cannot perform tissue response analysis.")
        return None
    
    if verbose:
        print("in prepare_tissue_response_kernel: duration_ms", duration_ms)

    std_dev_ms = plot_pms['tissue_response_kernel']['std_dev_ms']
    if verbose:
        print("in prepare_tissue_response_kernel: std_dev_ms", std_dev_ms)
    
    # Create the kernel
    tissue_response_kernel = gaussian_kernel(duration_ms, sampling_rate_Hz, std_dev_ms)
    if verbose:
        print("in prepare_tissue_response_kernel: len of kernel", len(tissue_response_kernel))
    if plot_pms['tissue_response_kernel']['plot']:
        kernel_plot(tissue_response_kernel, sampling_rate_Hz, plot_pms, 'tissue-like_response_')

    return tissue_response_kernel

def produce_exc_firing_rates_from_spike_times(nest_pms, plot_pms, crop_pms, analysis_pms, cropped_events, verbose):
    if verbose:
        print("in produce_exc_firing_rates_from_spike_times")
    # Collect spike times from all populations into a single vector 
    num_exc_pop = nest_pms['network']['num_exc_pop']
    single_trace_spike_times = combine_spike_times_in_single_trace(cropped_events, num_exc_pop)
    # Calculate combined firing rates
    num_exc_neu_per_pop = nest_pms['network']['num_exc_neu_per_pop']
    time_points, combined_firing_rate = calculate_firing_rate(
        crop_pms, analysis_pms, single_trace_spike_times, num_exc_pop * num_exc_neu_per_pop)
    if plot_pms['instantaneous_firing_rate']:
        # Plot of combined instantaneous firing rates   
        firing_rates_plot(time_points, combined_firing_rate, crop_pms, plot_pms, 'instantaneous_')
    return time_points, combined_firing_rate

def estimate_firing_rates_from_spike_like_waveforms(
    analysis_pms, plot_pms, crop_pms,
    time_points, combined_firing_rate, verbose):
    # Creating the kernel to smooth an instantaneous spike time into a spike-like waveform
    single_spike_kernel = prepare_single_spike_kernel(analysis_pms, plot_pms, crop_pms, verbose)
    if single_spike_kernel is None:
        print("Skipping spike-like waveform analysis due to short cropping duration.")
        return time_points, None, None
    # Prepare firing rate-like signals from sums of spike-like waveforms
    smoothed_spikes_firing_rate = smooth_signal(combined_firing_rate, single_spike_kernel)
    if verbose:
        print("in estimate_firing_rates_from_spike_like_waveforms: len of smoothed signal", len(smoothed_spikes_firing_rate))
        print("in estimate_firing_rates_from_spike_like_waveforms: type of smoothed signal", type(smoothed_spikes_firing_rate))
    # Producing the estimation of firing rate from spike-like waves 
    if plot_pms['smoothed_spikes_firing_rate']:
        firing_rates_plot(time_points, smoothed_spikes_firing_rate, crop_pms, plot_pms, 'spike-like_waveforms_') 
    return time_points, smoothed_spikes_firing_rate, single_spike_kernel

def produce_high_freq_spectrum_from_spike_like_waveforms(
    plot_pms, analysis_pms, smoothed_spikes_firing_rate, verbose):
    if smoothed_spikes_firing_rate is None:
        print("Cannot produce high-frequency spectrum: smoothed spikes firing rate is None.")
        return
    if verbose:
        print("in produce_high_freq_spectrum_from_spike_like_waveforms:")
    # Dividing the data into segments for reliable spectral analysis
    max_plot_freq_Hz = plot_pms['high_freq_spectrum']['max_plot_freq_Hz']
    smoothing_length = plot_pms['high_freq_spectrum']['smoothing_length']  # Frequency samples
    if plot_pms['high_freq_spectrum']['plot']:
        compute_spectrum_with_error_bands(
            smoothed_spikes_firing_rate, analysis_pms,
            max_plot_freq_Hz, smoothing_length, plot_pms, 'spike-like_waveforms_')
    return

def produce_firing_rates_from_spikes(nest_pms, plot_pms, crop_pms, analysis_pms, cropped_events, verbose):
    if verbose: 
        print("in produce_firing_rates_from_spikes")
    # Produce firing rates from spike times
    time_points, combined_firing_rate = produce_exc_firing_rates_from_spike_times(
        nest_pms, plot_pms, crop_pms, analysis_pms, cropped_events, verbose)
    # Estimating the firing rates from the spike-like waveforms
    # They are obtained using a spike-waform-like filter
    time_points, smoothed_spikes_firing_rate, single_spike_kernel = estimate_firing_rates_from_spike_like_waveforms(
        analysis_pms, plot_pms, crop_pms,
        time_points, combined_firing_rate, verbose)
    return time_points, smoothed_spikes_firing_rate

def produce_high_freq_spectrogram_from_spike_like_waveforms(
    plot_pms, analysis_pms, time_points, smoothed_spikes_firing_rate):
    if smoothed_spikes_firing_rate is None:
        print("Cannot produce high-frequency spectrogram: smoothed spikes firing rate is None.")
        return
    max_plot_freq_Hz = plot_pms['high_freq_spectrogram']['max_plot_freq_Hz']
    # Spectrogram plot
    if plot_pms['high_freq_spectrogram']['plot']:
        plot_spectrogram(
            time_points[0],
            smoothed_spikes_firing_rate,
            analysis_pms, max_plot_freq_Hz, plot_pms, 'spike-like_waveforms_')
    return

# This produces activity rates, spectrum, and spectrogram
def high_frequency_analysis_from_spikes(
    nest_pms, plot_pms, crop_pms, analysis_pms,
    cropped_events, verbose):
    # 1- ACTIVITY RATES
    # This includes a filter that transforms the spike times into spike-like waveforms
    time_points, smoothed_spikes_firing_rate = produce_firing_rates_from_spikes(
        nest_pms, plot_pms, crop_pms, analysis_pms, cropped_events, verbose)
    if smoothed_spikes_firing_rate is None:
        print("High-frequency analysis cannot be performed due to earlier issues.")
        return time_points, None
    # 2- SPECTRUM,
    # High-frequency spectrum with error bands
    produce_high_freq_spectrum_from_spike_like_waveforms(
        plot_pms, analysis_pms, smoothed_spikes_firing_rate, verbose)
    # 3- SPECTROGRAM
    # High-frequency spectrogram from spike-like waveforms
    produce_high_freq_spectrogram_from_spike_like_waveforms(
        plot_pms, analysis_pms, time_points, smoothed_spikes_firing_rate)
    return time_points, smoothed_spikes_firing_rate

def produce_tissue_response_activity_spectrum_spectrogram(
    nest_pms, plot_pms, crop_pms, analysis_pms,
    time_points, smoothed_spikes_firing_rate, verbose):
    if smoothed_spikes_firing_rate is None:
        print("Cannot perform tissue response analysis: smoothed spikes firing rate is None.")
        return
    # 1- KERNEL: creating a kernel to smooth from the spike-like waveforms 
    # towards tissue- and synaptic-like responses for ECG/LFP/ECoG-like analysis
    tissue_response_kernel = prepare_tissue_response_kernel(analysis_pms, plot_pms, crop_pms, verbose)
    if tissue_response_kernel is None:
        print("Skipping tissue response analysis due to short cropping duration.")
        return  # Exit the function early

    # 2- ACTIVITY evaluation
    # Evaluating activity from tissue response
    tissue_response_rate = smooth_signal(smoothed_spikes_firing_rate, tissue_response_kernel)
    if verbose: 
        print("len of tissue signal", len(tissue_response_rate))
        print("type of tissue signal", type(tissue_response_rate))
    if plot_pms['tissue_response_rate']:
        firing_rates_plot(time_points, tissue_response_rate, crop_pms, plot_pms, 'tissue-like_response_') 

    # 3- SPECTRUM
    # Dividing the data into segments for reliable spectral analysis
    max_plot_freq_Hz = plot_pms['tissue_response_spectrum']['max_plot_freq_Hz']
    smoothing_length = plot_pms['tissue_response_spectrum']['smoothing_length']  # Frequency samples
    if plot_pms['tissue_response_spectrum']['plot']:
        compute_spectrum_with_error_bands(
            tissue_response_rate, analysis_pms, max_plot_freq_Hz, smoothing_length, plot_pms, 'tissue-like_')
        max_plot_freq_Hz = plot_pms['tissue_response_spectrogram']['max_plot_freq_Hz']

    # 4- SPECTROGRAM
    if plot_pms['tissue_response_spectrogram']['plot']:
        plot_spectrogram(
            time_points[0],
            tissue_response_rate,
            analysis_pms, max_plot_freq_Hz, plot_pms, 'tissue-like_')
    return

def produce_rastegrams_rates_spectra_spectrograms(
    directories_and_list_of_yamls,
    nest_pms, crop_pms, plot_pms, analysis_pms,
    cropped_events, cropped_inh_events, verbose):
    # 0- Add to plot_pms info about the absolute path to possible storage of plots

    # 1- RASTEGRAMS
    # Produce the rastegrams
    # PLEASE: select if the exc, inh or both are needed from yaml files
    produce_rastegrams(nest_pms, plot_pms, crop_pms, cropped_events, cropped_inh_events, verbose)
    if plot_pms['intermediate_plt_show']:
        plt.show()
        print("FROM SPIKE-LIKE WAVEFORMS:")
        print("--------------------------")
    # 2- SPIKE-like waveforms to RATES, SPECTRUM and SPECTROGRAM 
    # Includes filter to transform spike-times into spike-like waveforms
    # PLEASE: use the yaml files to select the plots and the filter definition 
    time_points, smoothed_spikes_firing_rate = high_frequency_analysis_from_spikes(
        nest_pms, plot_pms, crop_pms, analysis_pms,
        cropped_events, verbose)
    if plot_pms['intermediate_plt_show']:
        plt.show()
        print("--------------------------")    
        print("FROM TISSUE-LIKE WAVEFORMS ")
        print("--------------------------")
    # 3- TISSUE-like to RATES, SPECTRUM and SPECTROGRAM  
    # Includes a kernel to smooth from the spike-like waveforms 
    # towards tissue- and synaptic-like responses for ECG/LFP/ECoG-like analysis
    # PLEASE: use the yaml file to decide the plots and kernel
    if smoothed_spikes_firing_rate is not None:
        produce_tissue_response_activity_spectrum_spectrogram(
            nest_pms, plot_pms, crop_pms, analysis_pms,
            time_points, smoothed_spikes_firing_rate, verbose)
    else:
        print("Skipping tissue response analysis due to previous issues.")
    plt.show()
    return


