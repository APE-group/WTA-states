#!/usr/bin/env python
# coding: utf-8

# Imports
import nest
import numpy as np
import matplotlib.pyplot as plt
from scipy.signal import convolve, spectrogram, butter, filtfilt

# In debug_mode, giving a first look before starting the analysis
def preliminary_sim_look(debug_mode, nest_pms, spike_recorders, inh_spike_recorder, recording_pms):
    if debug_mode:
        num_exc_neu_per_pop = nest_pms["network"]["num_exc_neu_per_pop"]
        num_exc_pop = nest_pms["network"]["num_exc_pop"]
        num_inh_neu = nest_pms["network"]["num_inh_neu"]
        total_exc_pop_spikes = [0 for _ in range(num_exc_pop)] 
        average_exc_pop_firing_rate_Hz = [0 for _ in range(num_exc_pop)]

        for i in range(num_exc_pop):
            d = nest.GetStatus(spike_recorders[i], "events")[0]
            times = d.get("times", [])
            senders = d.get("senders", [])
            total_exc_pop_spikes[i] = len(times)

            if total_exc_pop_spikes[i] > 0:
                print(f"pop {i} first recorded event at time {times[0]} from sender {senders[0]}")
                print(f"pop {i} last recorded event at time {times[-1]} from sender {senders[-1]}")
            else:
                print(f"pop {i} No spikes recorded.")

            print(f"pop {i} total_exc_pop_spikes = {total_exc_pop_spikes[i]}")
            average_exc_pop_firing_rate_Hz[i] = (
                total_exc_pop_spikes[i] / num_exc_neu_per_pop
            ) * (1000.0 / recording_pms["duration_ms"])  # 1000 ms in one s
            print(f"pop {i} average_exc_pop_firing_rate_Hz = {average_exc_pop_firing_rate_Hz[i]} Hz")

        d_inh = nest.GetStatus(inh_spike_recorder, "events")[0]
        times_inh = d_inh.get("times", [])
        senders_inh = d_inh.get("senders", [])
        total_inh_spikes = len(times_inh)

        if total_inh_spikes > 0:
            print(f"INHIBITORIES first recorded event at time {times_inh[0]} from sender {senders_inh[0]}")
            print(f"INHIBITORIES last recorded event at time {times_inh[-1]} from sender {senders_inh[-1]}")
        else:
            print("INHIBITORIES No spikes recorded.")

        print(f"INHIBITORIES total_inh_pop_spikes = {total_inh_spikes}")
        average_inh_firing_rate_Hz = (
            total_inh_spikes / num_inh_neu
        ) * (1000.0 / recording_pms["duration_ms"])  # 1000 ms in one s
        print(f"INHIBITORIES average_inh_firing_rate_Hz = {average_inh_firing_rate_Hz} Hz")

def crop_events_from_spike_recorders(crop_pms, spike_recorders):
    cropped_events = []  # Initialize an empty list to hold the cropped events for each recorder
    for recorder in spike_recorders:
        # Get events for the current recorder
        recorder_events = nest.GetStatus(recorder, "events")[0]
        spike_times = recorder_events.get("times", [])
        spike_senders = recorder_events.get("senders", [])
        
        # Filter events based on the specified crop start and stop times
        cropped_indices = [i for i, time in enumerate(spike_times) if crop_pms["start_ms"] <= time <= crop_pms["stop_ms"]]
        cropped_times = [spike_times[i] for i in cropped_indices]
        cropped_senders = [spike_senders[i] for i in cropped_indices]
        
        # Append the cropped events for the current recorder to the list
        cropped_events.append({
            "times": cropped_times,
            "senders": cropped_senders
        })
      
    return cropped_events

def crop_inh_events(crop_pms, inh_spike_recorder):
    cropped_events = []  # Initialize an empty list to hold the cropped events
    recorder_events = nest.GetStatus(inh_spike_recorder, "events")[0]
    spike_times = recorder_events.get("times", [])
    spike_senders = recorder_events.get("senders", [])
        
    # Filter events based on the specified crop start and stop times
    cropped_indices = [i for i, time in enumerate(spike_times) if crop_pms["start_ms"] <= time <= crop_pms["stop_ms"]]
    cropped_times = [spike_times[i] for i in cropped_indices]
    cropped_senders = [spike_senders[i] for i in cropped_indices]
        
    # Append the cropped events to the list
    cropped_events.append({
        "times": cropped_times,
        "senders": cropped_senders
    })
    return cropped_events

def derive_analysis_pms_from_sampling_pms(sampling_pms):
    analysis_pms = {}
    analysis_pms["spectrogram_window_ms"] = sampling_pms["spectrogram_window_ms"]
    analysis_pms["lower_detectable_frequency_Hz"] = (
        1 / (analysis_pms["spectrogram_window_ms"] / 1000.0)
    )  # 1000 from ms to s    
    analysis_pms["spikes_sampling_window_ms"] = sampling_pms["spikes_sampling_window_ms"]
    analysis_pms["analysis_sampling_freq_Hz"] = (
        1 / (analysis_pms["spikes_sampling_window_ms"] / 1000.0)
    )  # 1000 from ms to s
    analysis_pms["max_nyquist_Hz"] = analysis_pms["analysis_sampling_freq_Hz"] / 2
    analysis_pms["low_freq_sampling_window_ms"] = sampling_pms["low_freq_sampling_window_ms"]
    analysis_pms["low_freq_sampling_Hz"] = (
        1 / (analysis_pms["low_freq_sampling_window_ms"] / 1000.0)
    )  # 1000 from ms to s
    analysis_pms["max_low_freq_nyquist_Hz"] = analysis_pms["low_freq_sampling_Hz"] / 2    
    return analysis_pms

def analysis_pms_print(recording_pms, crop_pms, analysis_pms):
    print("---RECORDING INTERVAL---")
    print("start recording at time", recording_pms["start_ms"], "ms")
    print("stop recording at time", recording_pms["stop_ms"], "ms")
    print("recording_duration", recording_pms["duration_ms"], "ms")
    print("---CROPPING INTERVAL---")
    print("crop start at time", crop_pms["start_ms"], "ms")
    print("crop stop at time", crop_pms["stop_ms"], "ms")
    print("crop duration", crop_pms["duration_ms"], "ms")
    print("---HIGH FREQ analysis AND (spikes)---")
    print("spectrogram window", analysis_pms["spectrogram_window_ms"], "ms")
    print("lower_detectable_frequency", analysis_pms["lower_detectable_frequency_Hz"], "Hz")
    print("spikes sampling window", analysis_pms["spikes_sampling_window_ms"], "ms")
    print("analysis sampling freq", analysis_pms["analysis_sampling_freq_Hz"], "Hz")
    print("max nyquist frequency", analysis_pms["max_nyquist_Hz"], "Hz")
    print("----LOW FREQ analysis (ECG, LFP like)----")
    print("(low freq) spectrogram window", analysis_pms["spectrogram_window_ms"], "ms")
    print("lower_detectable_frequency", analysis_pms["lower_detectable_frequency_Hz"], "Hz")          
    print("low freq sampling window", analysis_pms["low_freq_sampling_window_ms"], "ms")
    print("low freq sampling freq", analysis_pms["low_freq_sampling_Hz"], "Hz")
    print("(low_freq) max nyquist frequency", analysis_pms["max_low_freq_nyquist_Hz"], "Hz")
    print("(low freq) spectrogram window", analysis_pms["spectrogram_window_ms"], "ms")

# Plot of rastergram of excitatories
def exc_pops_rastegram_plot(cropped_events, num_exc_pop, crop_pms, plot_pms):
    fig, ax = plt.subplots(figsize=(12, 8))
    colours = ["blue", "green", "orange", "black", "purple", "cyan", "blue", "green", "orange",
               "black", "purple", "cyan", "blue", "green", "orange", "black", "purple",
               "cyan", "blue", "green", "orange", "black", "purple", "cyan"]
    assert len(colours) >= num_exc_pop
    
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
    fig_file_name = plot_pms['absolute_plot_dir_path'] + '/' + 'exc_rastergram.pdf'
    if plot_pms['save_plots_in_output_dir']:
        plt.savefig(fig_file_name)
    # plt.show()

# Plot of rastergram of inhibitories
def inh_rastergram_plot(cropped_inh_events, crop_pms, plot_pms):
    fig, ax = plt.subplots(figsize=(12, 4))
    # Plot for inhibitory neurons
    inh_pop_events = cropped_inh_events[0]
    ax.scatter(inh_pop_events["times"], inh_pop_events["senders"], s=5, c="black", label="Inhibitory Population")
    ax.set_xlim(crop_pms["start_ms"], crop_pms["stop_ms"])
    ax.set_title("Rastergram of Inhibitory Neural Activity")
    ax.set_xlabel("Time (ms)")
    ax.set_ylabel("Neuron ID")
    ax.legend()
    plt.tight_layout()
    fig_file_name = plot_pms['absolute_plot_dir_path'] + '/' + 'inh_rastergram.pdf'
    if plot_pms['save_plots_in_output_dir']:
        plt.savefig(fig_file_name)
    # plt.show()

def combine_spike_times_in_single_trace(cropped_events, num_pops):
    combined_spike_times = []
    for pop in range(num_pops):
        pop_events = cropped_events[pop]
        combined_spike_times.extend(pop_events["times"])
    # Ensure the spike times are sorted as np.histogram expects ordered data for bin edges
    combined_spike_times.sort()
    return combined_spike_times                          

def calculate_firing_rate(crop_pms, analysis_pms, spike_times, num_neu):
    local_debug = False
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
        print("number of expected time intervals sampled", len(time_points) - 1)
    assert len(time_points) >= 2
    # Using time_points as bin edges for the histogram
    spike_counts, bins = np.histogram(spike_times, bins=time_points)
    if local_debug:
        print("length of downsampled vector", len(spike_counts))
    # Calculate firing rate: adjust for number of neurons and convert count to rate per second
    firing_rate = (spike_counts / num_neu) * (1000.0 / analysis_pms["spikes_sampling_window_ms"])

    # Since we use time_points as bins edges, the time for each firing rate value should be the midpoint of each bin
    firing_time_points = (bins[:-1] + bins[1:]) / 2

    return firing_time_points, firing_rate


def firing_rates_plot(time_points, firing_rate, crop_pms, plot_pms, title_preamble_string):
    # Ensure lengths are consistent
    if len(time_points) != len(firing_rate):
        min_length = min(len(time_points), len(firing_rate))
        time_points = time_points[:min_length]
        firing_rate = firing_rate[:min_length]
    
    fig, ax = plt.subplots(figsize=(12, 4))
    ax.plot(time_points, firing_rate)
    ax.set_xlim(crop_pms["start_ms"], crop_pms["stop_ms"])
    ax.set_xlabel("Time (ms)")
    ax.set_ylabel("Firing Rate (Hz)")
    ax.set_title("Firing Rates")
    fig_file_name = plot_pms['absolute_plot_dir_path'] + '/' + title_preamble_string + 'firing_rates.pdf'
    if plot_pms['save_plots_in_output_dir']:
        plt.savefig(fig_file_name)
    # plt.show()


def gaussian_kernel(duration_ms, sampling_rate_Hz, std_dev_ms):
    local_debug = False
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
    sampling_duration_ms = sampling_duration_s * 1000.0  # ms / s
    kernel_start_ms = -duration_ms / 2.0
    kernel_stop_ms = duration_ms / 2.0
    if local_debug:
        print("kernel_start_ms", kernel_start_ms)
        print("kernel_stop_ms", kernel_stop_ms)
        print("sampling_duration_ms", sampling_duration_ms)       
    t_ms = np.arange(kernel_start_ms, kernel_stop_ms, sampling_duration_ms)
    if local_debug:
        print("length in samples", len(t_ms))
    gaussian = np.exp(-0.5 * (t_ms / std_dev_ms) ** 2)
    gaussian /= gaussian.sum()  # Normalize the kernel to ensure the signal energy is preserved
    return gaussian

# Plot of kernel
def kernel_plot(kernel, sampling_rate_Hz, plot_pms, title_preamble_string):
    sampling_duration_s = 1 / sampling_rate_Hz
    sampling_duration_ms = sampling_duration_s * 1000.0  # ms / s
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
    fig_file_name = plot_pms['absolute_plot_dir_path'] + '/' + title_preamble_string + 'kernel.pdf'
    if plot_pms['save_plots_in_output_dir']:
        plt.savefig(fig_file_name)
    # plt.show()

def smooth_signal(signal, kernel):
    """
    Smooths the given signal by convolving it with the specified kernel.
    
    Parameters:
    - signal: The original signal 
    - kernel: The smoothing kernel.
    
    Returns:
    - Smoothed signal.
    """
    return convolve(signal, kernel, mode='same')

def moving_average(data, window_size):
    window = np.ones(int(window_size)) / float(window_size)
    smoothed = np.convolve(data, window, 'same')
    edge_correction = np.convolve(np.ones_like(data), window, 'same')
    return smoothed / edge_correction

def compute_spectrum_segments(data, sampling_rate, segment_length_ms, overlap_factor=0.5):
    segment_length = int(segment_length_ms * sampling_rate / 1000)
    step = int(segment_length * (1 - overlap_factor))
    num_segments = (len(data) - segment_length) // step + 1

    if num_segments <= 0:
        print("Warning: Not enough data segments for spectrum computation.")
        return np.array([]), np.array([]), 0

    spectra = []

    for i in range(num_segments):
        segment = data[i * step: i * step + segment_length]
        if len(segment) < segment_length:
            # Pad the segment if it's shorter than expected
            segment = np.pad(segment, (0, segment_length - len(segment)), 'constant')
        fft_data = np.fft.fft(segment)
        fft_freq = np.fft.fftfreq(segment_length, d=1 / sampling_rate)
        magnitude = np.abs(fft_data)
        spectra.append(magnitude)

    return np.array(spectra), fft_freq, num_segments

def plot_spectrum_with_error_bands(
    segment_duration_ms,
    max_plot_freq_Hz,
    spectra,
    freqs,
    num_segments,
    smoothing_length,
    plot_pms,
    title_preamble_string,
):
    mean_spectrum = np.mean(spectra, axis=0)
    smoothed_mean_spectrum = moving_average(mean_spectrum, smoothing_length)

    # Ensure consistent array lengths
    min_length = min(len(freqs), len(mean_spectrum), len(smoothed_mean_spectrum))
    freqs = freqs[:min_length]
    mean_spectrum = mean_spectrum[:min_length]
    smoothed_mean_spectrum = smoothed_mean_spectrum[:min_length]

    assert segment_duration_ms > 0.0, "Segment duration must be positive."
    lower_detectable_frequency_Hz = 1000.0 / segment_duration_ms

    valid_indices = (freqs >= lower_detectable_frequency_Hz) & (freqs <= max_plot_freq_Hz)
    freqs = freqs[valid_indices]
    mean_spectrum = mean_spectrum[valid_indices]
    smoothed_mean_spectrum = smoothed_mean_spectrum[valid_indices]

    if len(freqs) == 0:
        print("Warning: No valid frequency data after applying filters.")
        return

    max_index = np.argmax(mean_spectrum)  # Index of the max value in the mean spectrum
    max_freq_Hz = freqs[max_index]  # Corresponding frequency

    plt.figure(figsize=(12, 8))
    # Use try-except to handle cases where percentile computation fails
    try:
        lower_quartile = np.percentile(spectra, 25, axis=0)[valid_indices]
        upper_quartile = np.percentile(spectra, 75, axis=0)[valid_indices]
        plt.fill_between(freqs, lower_quartile, upper_quartile, color='lightblue', label='Quartile Range')
    except IndexError:
        print("Warning: Unable to compute quartiles for spectra.")
        lower_quartile = upper_quartile = np.zeros_like(mean_spectrum)

    plt.plot(freqs, mean_spectrum, label='Mean Spectrum', color='blue')
    plt.plot(freqs, smoothed_mean_spectrum, label='Smoothed Mean Spectrum', color='red', linestyle='--')
    plt.semilogy()
    plt.xlim(lower_detectable_frequency_Hz, max_plot_freq_Hz)
    plt.xlabel('Frequency (Hz)')
    plt.ylabel('Magnitude (Log Scale)')
    # Frequency resolution in Hz
    smoothing_freq_Hz = lower_detectable_frequency_Hz * smoothing_length
    plt.title(f'Spectrum with Error Bands\nNumber of segments: {num_segments}, Segment duration: {segment_duration_ms} ms, Smoothing Window: {smoothing_freq_Hz:.2f} Hz')
    plt.legend()
    plt.grid(True)

    # Annotating the max frequency at the bottom of the plot
    plt.annotate(
        f'Max freq: {max_freq_Hz:.2f} Hz',
        xy=(max_freq_Hz, smoothed_mean_spectrum[max_index]),
        xytext=(max_freq_Hz, 0.1 * np.min(smoothed_mean_spectrum)),
        arrowprops=dict(facecolor='black', shrink=0.05, width=1, headwidth=8),
        verticalalignment='bottom',
    )
    fig_file_name = plot_pms['absolute_plot_dir_path'] + '/' + title_preamble_string + 'spectrum_with_error_bands.pdf'
    if plot_pms['save_plots_in_output_dir']:
        plt.savefig(fig_file_name)

    # plt.show()

def compute_spectrum_with_error_bands(data, analysis_pms, max_plot_freq_Hz, smoothing_length, plot_pms, title_preamble_string):
    segment_duration_ms = analysis_pms["spectrogram_window_ms"]
    analysis_sampling_freq_Hz = analysis_pms["analysis_sampling_freq_Hz"]

    # Check if data is sufficient for spectrum computation
    if len(data) < 2:
        print(f"Warning: Not enough data points ({len(data)}) for spectrum analysis.")
        return

    # Compute spectra
    spectra, fft_freq, num_segments = compute_spectrum_segments(
        data, analysis_sampling_freq_Hz, segment_duration_ms
    )

    # Check if spectra were computed successfully
    if spectra.size == 0 or fft_freq.size == 0:
        print("Warning: Spectra computation resulted in empty arrays.")
        return

    # Ensure that smoothed_mean_spectrum has the same length as mean_spectrum
    mean_spectrum = np.mean(spectra, axis=0)
    smoothed_mean_spectrum = moving_average(mean_spectrum, smoothing_length)

    # Truncate smoothed_mean_spectrum to match mean_spectrum length if necessary
    min_length = min(len(mean_spectrum), len(smoothed_mean_spectrum))
    mean_spectrum = mean_spectrum[:min_length]
    smoothed_mean_spectrum = smoothed_mean_spectrum[:min_length]
    fft_freq = fft_freq[:min_length]

    plot_spectrum_with_error_bands(
        segment_duration_ms,
        max_plot_freq_Hz,
        spectra,
        fft_freq,
        num_segments,
        smoothing_length,
        plot_pms,
        title_preamble_string,
    )

# Spectrogram plot functions
def plot_spectrogram(time_origin_ms, data, analysis_pms, max_plot_freq_Hz, plot_pms, title_preamble_string):
    lower_detectable_freq_Hz = analysis_pms["lower_detectable_frequency_Hz"]
    max_nyquist_Hz = analysis_pms["max_nyquist_Hz"]
    if max_plot_freq_Hz > max_nyquist_Hz:
        print(f"Warning: max_plot_freq_Hz ({max_plot_freq_Hz} Hz) exceeds Nyquist frequency ({max_nyquist_Hz} Hz). Adjusting max_plot_freq_Hz.")
        max_plot_freq_Hz = max_nyquist_Hz
    fig, ax = plt.subplots(figsize=(12, 4))
    # Compute the number of data points per segment based on the sampling frequency and window size in ms
    nperseg = int(analysis_pms["spectrogram_window_ms"] * analysis_pms["analysis_sampling_freq_Hz"] / 1000.0)

    # Calculate the spectrogram
    f, t, Sxx = spectrogram(data, fs=analysis_pms["analysis_sampling_freq_Hz"], nperseg=nperseg)
    scaled_shifted_t = t * 1000 + time_origin_ms  # Bringing to ms
    
    # Use actual time_points for the x-axis
    c = ax.pcolormesh(scaled_shifted_t, f, 10 * np.log10(Sxx + 1e-10), shading='gouraud')  # +1e-10 to avoid log(0)
    ax.set_ylim(lower_detectable_freq_Hz, max_plot_freq_Hz)
    ax.set_ylabel('Frequency [Hz]')
    ax.set_xlabel('Time [ms]')
    plt.title('Spectrogram')

    # Add a color bar on the right to show the numerical values of different colors
    plt.colorbar(c, ax=ax, label='Intensity (dB)')
    fig_file_name = plot_pms['absolute_plot_dir_path'] + '/' + title_preamble_string + 'spectrogram.pdf'
    if plot_pms['save_plots_in_output_dir']:
        plt.savefig(fig_file_name)

    # plt.show()

def set_final_crop_pms(crop_pms, nest_pms, verbose):
    # Setting the cropping parameters
    default_cropping = crop_pms['default_cropping']
    
    if default_cropping:
        crop_pms['start_ms'] = nest_pms["recording_pms"]["start_ms"]
        crop_pms['stop_ms'] = nest_pms["recording_pms"]["stop_ms"]
        crop_pms['duration_ms'] = crop_pms["stop_ms"] - crop_pms["start_ms"]
        crop_pms['delay_ms'] = 0.0    
    else:
        crop_pms['start_ms'] = nest_pms["recording_pms"]["start_ms"] + crop_pms['delay_ms']
        assert(crop_pms['delay_ms']>=0.0 and crop_pms['start_ms']<=nest_pms["recording_pms"]["stop_ms"])
        crop_pms['stop_ms'] = crop_pms['start_ms'] + crop_pms['duration_ms']
        assert(crop_pms['duration_ms']>=0.0 and crop_pms['stop_ms']<=nest_pms["recording_pms"]["stop_ms"])
        crop_pms['duration_ms'] = crop_pms["stop_ms"] - crop_pms["start_ms"]
        
    
    if verbose:
        print("FINAL crop_pms values", crop_pms)
    return crop_pms

# Decision of spectral_window_ms based on the value of the crop_duration_ms
def spectral_window_ms_analysis_switch(crop_duration_ms):
    if crop_duration_ms < 4000.0:
        spectral_window_ms = crop_duration_ms / 4.0
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

# Function to produce firing rates from spikes (included in the module)
def produce_firing_rates_from_spikes(nest_pms, plot_pms, crop_pms, analysis_pms, cropped_events, verbose):
    num_exc_neu_per_pop = nest_pms["network"]["num_exc_neu_per_pop"]
    num_exc_pop = nest_pms["network"]["num_exc_pop"]

    # Combine spike times from all excitatory populations
    combined_spike_times = combine_spike_times_in_single_trace(cropped_events, num_exc_pop)
    if len(combined_spike_times) == 0:
        print("Warning: No spikes recorded in excitatory populations.")
        return np.array([]), np.array([])

    # Calculate firing rate using the total number of neurons
    time_points, firing_rate = calculate_firing_rate(
        crop_pms, analysis_pms, combined_spike_times, num_exc_neu_per_pop * num_exc_pop
    )

    # Plot the firing rates
    if plot_pms.get('firing_rates_plot', {}).get('plot', False):
        firing_rates_plot(time_points, firing_rate, crop_pms, plot_pms, 'spike-like_waveforms_')

    # Generate and plot the kernel
    if plot_pms.get('kernel', {}).get('plot', False):
        kernel = gaussian_kernel(
            plot_pms['kernel']['duration_ms'],
            analysis_pms["analysis_sampling_freq_Hz"],
            plot_pms['kernel']['std_dev_ms']
        )
        kernel_plot(kernel, analysis_pms["analysis_sampling_freq_Hz"], plot_pms, 'spike-like_waveforms_')
    else:
        kernel = gaussian_kernel(
            plot_pms['kernel']['duration_ms'],
            analysis_pms["analysis_sampling_freq_Hz"],
            plot_pms['kernel']['std_dev_ms']
        )

    # Smooth the firing rate signal
    smoothed_firing_rate = smooth_signal(firing_rate, kernel)
    return time_points, smoothed_firing_rate