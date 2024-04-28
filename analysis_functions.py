#!/usr/bin/env python
# coding: utf-8

#imports
import numpy as np
import matplotlib.pyplot as plt
from scipy.signal import spectrogram, butter, filtfilt
    
def derive_analysis_pms_from_sampling_pms(sampling_pms):
    analysis_pms={}
    analysis_pms["spectrogram_window_ms"]=sampling_pms["spectrogram_window_ms"]
    analysis_pms["lower_detectable_frequency_Hz"] =\
        1 / (analysis_pms["spectrogram_window_ms"] / 1000.0) # 1000 from ms to s    
    analysis_pms["spikes_sampling_window_ms"]=sampling_pms["spikes_sampling_window_ms"] 
    analysis_pms["analysis_sampling_freq_Hz"] =\
        1 / (analysis_pms["spikes_sampling_window_ms"] /1000.0) # 1000 from ms to s
    analysis_pms["max_nyquist_Hz"] = analysis_pms["analysis_sampling_freq_Hz"]/2
    analysis_pms["low_freq_sampling_window_ms"] = sampling_pms["low_freq_sampling_window_ms"]
    analysis_pms["low_freq_sampling_Hz"]=\
        1 / (analysis_pms["low_freq_sampling_window_ms"] /1000.0) # 1000 from ms to s
    analysis_pms["max_low_freq_nyquist_Hz"] = analysis_pms["low_freq_sampling_Hz"]/2    
    return analysis_pms

def analysis_pms_print(recording_pms,crop_pms,analysis_pms):
    print("---RECORDING INTERVAL---")
    print("start recording at time", recording_pms["start_ms"], "ms")
    print("stop recording at time", recording_pms["stop_ms"], "ms")
    print("recording_duration", recording_pms["duration_ms"] ,"ms")
    print("---CROPPING INTERVAL")
    print("crop start at time", crop_pms["start_ms"], "ms")
    print("crop stop at time", crop_pms["stop_ms"], "ms")
    print("crop duration",crop_pms["duration_ms"] ,"ms")
    print("---HIGH FREQ analysis AND (spikes) ")
    print("lower_detectable_frequency", analysis_pms["lower_detectable_frequency_Hz"], "Hz")
    print("spikes sampling window", analysis_pms["spikes_sampling_window_ms"], "ms")
    print("analysis sampling freq", analysis_pms["analysis_sampling_freq_Hz"], "Hz")
    print("max nyquist frequency", analysis_pms["max_nyquist_Hz"], "Hz")
    print("spectogram window", analysis_pms["spectrogram_window_ms"], "ms")
    print("----LOW FREQ analysis (ECG, LFP like)")
    print("lower_detectable_frequency", analysis_pms["lower_detectable_frequency_Hz"], "Hz")          
    print("low freq sampling window", analysis_pms["low_freq_sampling_window_ms"], "ms")
    print("low freq sampling freq", analysis_pms["low_freq_sampling_Hz"], "Hz")
    print("(low_freq) max nyquist frequency", analysis_pms["max_low_freq_nyquist_Hz"], "Hz")
    print("(low freq) spectrogram window", analysis_pms["spectrogram_window_ms"], "ms")

# plot of rastegram of excitatories
def exc_pops_rastegram_plot(cropped_events,num_exc_pop, crop_pms):
    fig, ax = plt.subplots(figsize=(12, 4))
    colours = ["blue","green","orange","black","purple","cyan","blue","green","orange",
               "black","purple","cyan","blue","green","orange","black","purple",
               "cyan","blue","green","orange","black","purple","cyan"]
    assert(len(colours)>=num_exc_pop)
    
    # Plot for excitatory neurons
    for pop in range(num_exc_pop):
        pop_events = cropped_events[pop]
        ax.scatter(pop_events["times"], pop_events["senders"], s=5, c=colours[pop], label=f"Population {pop + 1}")
    ax.set_xlim(crop_pms["start_ms"], crop_pms["stop_ms"])
    ax.set_title("Rastergram of Excitatory Neural Activity")
    ax.set_xlabel("Time (ms)")
    ax.set_ylabel("Neuron ID")
    ax.legend()
    plt.tight_layout()
    plt.show()

# plot of rastegram of inhibitories
def inh_rastergram_plot(cropped_inh_events, crop_pms):
    fig, ax = plt.subplots(figsize=(12, 4))
    # Plot for inhibitory neurons (there is a single population of inh, but it is necessary to extact it (sic!))   
    inh_pop_events=cropped_inh_events[0]
    ax.scatter(inh_pop_events["times"], inh_pop_events["senders"], s=5, c="black", label="Inhibitory Population")
    ax.set_xlim(crop_pms["start_ms"], crop_pms["stop_ms"])
    ax.set_title("Rastergram of Inhibitory Neural Activity")
    ax.set_xlabel("Time (ms)")
    ax.set_ylabel("Neuron ID")
    ax.legend()
    plt.tight_layout()
    plt.show()
    
def combine_spike_times_in_single_trace(cropped_events, num_pops):
    combined_spike_times = []
    for pop in range(num_pops):
        pop_events = cropped_events[pop]
        combined_spike_times.extend(pop_events["times"])
    # Ensure the spike times are sorted as np.histogram expects ordered data for bin edges
    combined_spike_times.sort()
    return combined_spike_times                          
                              
def calculate_firing_rate(crop_pms, analysis_pms, spike_times, num_neu):
    local_debug=False
    start_ms = crop_pms["start_ms"]
    stop_ms = crop_pms["stop_ms"]
    spikes_sampling_window_ms = analysis_pms["spikes_sampling_window_ms"]        
    if local_debug:
        print("start_ms", start_ms)
        print("stop_ms", stop_ms)
        print("spikes_sampling_window_ms", spikes_sampling_window_ms)
    
    # Define the edges of bins based on the specified sampling window
    time_points = np.arange(start_ms, stop_ms, spikes_sampling_window_ms)
    if local_debug:
        print("number of expected time intervals sampled", len(time_points)-1)
    assert(len(time_points)>=2)
    # Using time_points as bin edges for the histogram
    spike_counts, bins = np.histogram(spike_times, bins=time_points)
    if local_debug:
        print("length of downsampled vector", len(spike_counts))
    # Calculate firing rate: adjust for number of neurons and convert count to rate per second
    firing_rate = (spike_counts / num_neu) \
        * (1000.0 / analysis_pms["spikes_sampling_window_ms"])

    # Since we use time_points as bins edges, the time for each firing rate value should be the midpoint of each bin for better accuracy
    firing_time_points = (bins[:-1] + bins[1:]) / 2

    return firing_time_points, firing_rate

#plot combined firing rates
def firing_rates_plot(time_points, firing_rate, crop_pms):
    fig, ax = plt.subplots(figsize=(4, 2))
    ax.plot(time_points, firing_rate)
    ax.set_xlim(crop_pms["start_ms"], crop_pms["stop_ms"])
    ax.set_xlabel("Time (ms)")
    ax.set_ylabel("Firing Rate (Hz)")
    ax.set_title("Firing Rates")
    plt.show()
    
def gaussian_kernel(duration_ms, sampling_rate_Hz, std_dev_ms):
    local_debug=False
    """
    Generates a Gaussian kernel for smoothing.
    
    Parameters:
    - duration_ms: The duration over which the kernel is non-zero.
    - sampling_rate_Hz: The number of samples per second.
    - std_dev_ms: The standard deviation of the Gaussian kernel in ms.
    
    Returns:
    - Gaussian kernel array.
    """
    sampling_duration_s = 1 / sampling_rate_Hz
    sampling_duration_ms = sampling_duration_s * 1000.0 # ms / s
    kernel_start_ms = -duration_ms / 2.0
    kernel_stop_ms = duration_ms / 2.0
    if local_debug:
        print("kernel_start_ms", kernel_start_ms)
        print("kernel_stop_ms", kernel_stop_ms)
        print("sampling_duration_ms", sampling_duration_ms)       
    t_ms = np.arange(kernel_start_ms, kernel_stop_ms, sampling_duration_ms)
    if local_debug:
        print("length in samples",len(t_ms))
    gaussian = np.exp(-0.5 * (t_ms / std_dev_ms) ** 2)
    gaussian /= gaussian.sum()  # Normalize the kernel to ensure the signal energy is preserved
    return gaussian

#plot of kernel
def kernel_plot(kernel, sampling_rate_Hz):
    sampling_duration_s = 1 / sampling_rate_Hz
    sampling_duration_ms = sampling_duration_s * 1000.0 # ms / s
    kernel_length_samples = len(kernel)
    duration_ms = kernel_length_samples * sampling_duration_ms
    start_time_ms = -duration_ms / 2.0
    stop_time_ms = duration_ms / 2.0
    t_ms = np.arange(start_time_ms, stop_time_ms, sampling_duration_ms)
    fig, ax = plt.subplots(figsize=(4, 2))
    ax.plot(t_ms, kernel)
    ax.set_xlim(start_time_ms, stop_time_ms)
    ax.set_xlabel("Time (ms)")
    ax.set_ylabel("Kernel amplitude (density / ms)")
    ax.set_title("Kernel")
    plt.show()
    
def smooth_signal(signal, kernel):
    """
    Smooths the given signal by convolving it with the specified kernel.
    
    Parameters:
    - signal: The original signal 
    - kernel: The smoothing kernel.
    
    Returns:
    - Smoothed signal.
    """
    return np.convolve(signal, kernel, mode='same')

#Spectrogram plot functions
def plot_spectrogram(time_origin_ms, data, analysis_pms, max_plot_freq_Hz):
    fig, ax = plt.subplots(figsize=(18, 8))
    # Compute the number of data points per segment based on the sampling frequency and window size in ms
    nperseg = int(analysis_pms["spectrogram_window_ms"] * analysis_pms["analysis_sampling_freq_Hz"] / 1000.0)

    # Calculate the spectrogram
    f, t, Sxx = spectrogram(data, fs=analysis_pms["analysis_sampling_freq_Hz"], nperseg=nperseg)
    scaled_shifted_t = t * 1000 + time_origin_ms # bringing to ms
    
    # Use actual time_points for the x-axis
    # Note: The pcolormesh function will automatically handle the alignment of these points with the spectrogram data
    c = ax.pcolormesh(scaled_shifted_t, f, 10 * np.log10(Sxx + 1), shading='gouraud')  # +1 to avoid log(0) issues
    ax.set_ylim(0,max_plot_freq_Hz)
    ax.set_ylabel('Frequency [Hz]')
    ax.set_xlabel('Time [ms]')
    plt.title('Spectrogram')

    # Add a color bar on the right to show the numerical values of different colors
    plt.colorbar(c, ax=ax, label='Intensity (dB)')

    plt.show()
    

def moving_average(data, window_size):
    """Compute the moving average using a simple convolution."""
    window = np.ones(int(window_size)) / float(window_size)
    return np.convolve(data, window, 'same')

def compute_and_plot_spectrum_with_moving_average(data, sampling_rate_Hz, freq_window_size_Hz, max_freq_plot_Hz):
    """
    Computes the Fourier Transform to find the spectrum of the given data,
    calculates a moving average of the spectrum, and plots both.
    
    Parameters:
    - data: 1D array of time series data.
    - sampling_rate: Sampling frequency of the data in Hz.
    - window_size: Window size for the moving average in number of samples.
    """
    # Compute the FFT (Fast Fourier Transform)
    n = len(data)
    fft_data = np.fft.fft(data)
    fft_freq = np.fft.fftfreq(n, d=1/sampling_rate_Hz)

    # Compute the magnitude of the FFT
    magnitude = np.abs(fft_data)
    positive_frequencies = fft_freq > 0
    
    # Calculate moving average of the magnitude spectrum
    moving_avg_magnitude = moving_average(magnitude, freq_window_size_Hz)

    # Plotting the spectrum and its moving average
    plt.figure(figsize=(10, 5))
    plt.xlim(0, max_freq_plot_Hz)
    plt.semilogy(fft_freq[positive_frequencies], magnitude[positive_frequencies], label='Original Spectrum')
    plt.semilogy(fft_freq[positive_frequencies], moving_avg_magnitude[positive_frequencies], label='Moving Average', linestyle='--')
    plt.title('Power Spectrum and Moving Average on Logarithmic Scale')
    plt.xlabel('Frequency (Hz)')
    plt.ylabel('Magnitude (Log Scale)')
    plt.legend()
    plt.grid(True)
    plt.show()
    max_index = np.argmax(magnitude[positive_frequencies])
    max_freq = fft_freq[positive_frequencies][max_index]
    print("max of fft at", max_freq)
    
def low_pass_filter(data, sampling_rate, cutoff_freq):
    """
    Applies a low-pass Butterworth filter to the given data.
    
    Parameters:
    - data: 1D array of time series data.
    - sampling_rate: Sampling frequency of the data in Hz.
    - cutoff_freq: Cutoff frequency in Hz. Frequencies higher than this will be attenuated.
    
    Returns:
    - filtered_data: 1D array of the filtered time series data.
    """
    # Normalize the frequency by the Nyquist frequency (half the sampling rate)
    nyquist = 0.5 * sampling_rate
    norm_cutoff_freq = cutoff_freq / nyquist

    # Design the Butterworth filter
    b, a = butter(N=5, Wn=norm_cutoff_freq, btype='low', analog=False)

    # Apply the filter to the data using filtfilt, which applies the filter forward and backward
    filtered_data = filtfilt(b, a, data)

    return filtered_data




