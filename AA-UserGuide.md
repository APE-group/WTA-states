# User Guide for WTA-states

## Disclaimer

This software is still in preliminary stages of development. Feel free to contact us if you like to contribute.

## Introduction

Welcome to WTA-states! This User Guide will help you get started and walk you through the key features of the software. WTA-states involves exercises combining brain states and a soft winner-takes-all mechanism, making it a valuable tool for neuroscientific research and modeling. Whether you're a beginner or an advanced user, you'll find everything you need here to make the most of your experience.

## Installation

### Prerequisites

- Operating System: Linux or macOS
- Main Dependencies: Python 3.8+, NEST Simulator, Jupyter Lab
- Other dependencies: yaml, numpy, matplotlib, scipy.signal 
- System Requirements: Minimum 8GB RAM, CPU with multiple cores recommended

## Key Features
- MAIN FILE TO BE LAUNCHED: WTA-awake-sleep.ipynb

### Feature 1: Brain State Simulation

- **Description**: This feature allows you to simulate different brain states using configurable parameters defined in YAML files.
- **Usage**:
  - Step 1: Load the appropriate YAML configuration file for the simulation.
  - Step 2: Run the simulation script to generate results.

### Feature 2: Winner-Takes-All Mechanism

- **Description**: This feature helps simulate a soft winner-takes-all mechanism in neural assemblies.
- **Usage**:
  - Step 1: Configure the synaptic connections using YAML files.
  - Step 2: Execute the main simulation to observe the winner-takes-all behavior.

### Feature 3: Plasticity and Neural Parameters

- **Description**: The software provides options to simulate synaptic plasticity and adjust neural parameters.
- **Usage**:
  - Step 1: Modify the plasticity settings in the configuration files.
  - Step 2: Run the simulation to observe changes in network dynamics.

### Feature 4: Compare the behaviour of single and two compartment neurons
- **Description**: The software provides options to compare single and two compartment neurons, in particular during learning of novel memories.
- **Usage**:
  - Step 1: select the appropriate excitatory neuron model
  - Step 2: set the time window for learning new memories

## YAML Configuration Files

WTA-states uses several YAML files to configure different aspects of the simulation. Below is an explanation of each type of YAML file, along with details on the key parameters that can be specified.

### 1. `basic_general_config.yaml`
- **Purpose**: This file contains general configuration settings for the simulation, such as global parameters that apply to all aspects of the network, including neuron models, synapse settings, and recording options.
- **Key Parameters**:
  - `brain_state`: Defines the brain state for the simulation (e.g., awake, NREM).
  - `use_single_compartment_environment`: Boolean flag to indicate if a single-compartment environment is used.
  - `inh_neu_params_filename`: Specifies the filename for inhibitory neuron parameters.
  - `exc_neu_params_filename`: Specifies the filename for excitatory neuron parameters.
  - `cf`: Configuration factor that adjusts simulation parameters.
  - `network`: Network-specific configurations.
    - `network.num_exc_neu_per_pop`: Number of excitatory neurons per population.
    - `network.num_exc_pop`: Number of excitatory populations.
    - `network.conn_rule`: Connection rule for the network.
    - `network.p_conn_exc`: Probability of connection for excitatory synapses.
    - `network.p_conn_inh`: Probability of connection for inhibitory synapses.
    - `network.use_exc_recurrency`: Boolean flag to enable or disable excitatory recurrency.
    - `network.use_poisson_generators`: Boolean flag to enable or disable Poisson generators.
    - `network.inter_pop_conn`: Specifies inter-population connectivity settings.
    - `network.num_poisson_generators`: Number of Poisson generators used.
    - `network.use_dc_exc_injectors`: Boolean flag to enable direct current injectors for excitatory neurons.
    - `network.use_dc_inh_injector`: Boolean flag to enable direct current injectors for inhibitory neurons.
    - `network.default_plasticity`: Default plasticity settings.
      - `network.default_plasticity.synapse_model`: Synapse model used for plasticity.
      - `network.default_plasticity.tau_plus`: Time constant for potentiation.
      - `network.default_plasticity.lambda`: Learning rate for plasticity.
      - `network.default_plasticity.alpha`: Scaling factor for synaptic changes.
      - `network.default_plasticity.mu_plus`: Positive modulation parameter.
      - `network.default_plasticity.mu_minus`: Negative modulation parameter.
      - `network.default_plasticity.weight`: Initial synaptic weight.
      - `network.default_plasticity.Wmax`: Maximum allowable weight.
      - `network.default_plasticity.delay`: Synaptic delay.
  - `contextual_poisson`: Contextual Poisson generator settings.
    - `contextual_poisson.awake`: Settings for awake state.
      - `contextual_poisson.awake.basic_rate`: Basic firing rate for awake state.
      - `contextual_poisson.awake.poisson_weight`: Weight of Poisson input.
      - `contextual_poisson.awake.spreading_factor`: Spreading factor for Poisson input.
      - `contextual_poisson.awake.events`: List of events during awake state.
    - `contextual_poisson.NREM`: Settings for NREM state.
      - `contextual_poisson.NREM.basic_rate`: Basic firing rate for NREM state.
      - `contextual_poisson.NREM.poisson_weight`: Weight of Poisson input.
      - `contextual_poisson.NREM.spreading_factor`: Spreading factor for Poisson input.
      - `contextual_poisson.NREM.events`: List of events during NREM state.
  - `poisson`: General Poisson input settings.
    - `poisson.awake`: Settings for awake state.
      - `poisson.awake.poisson_basic_rate`: Basic rate for Poisson input.
      - `poisson.awake.poisson_delta`: Delta parameter for Poisson input.
      - `poisson.awake.poisson_noise_spreading_factor`: Noise spreading factor.
      - `poisson.awake.poisson_weight`: Weight of Poisson input.
    - `poisson.NREM`: Settings for NREM state.
      - `poisson.NREM.poisson_basic_rate`: Basic rate for Poisson input.
      - `poisson.NREM.poisson_delta`: Delta parameter for Poisson input.
      - `poisson.NREM.poisson_noise_spreading_factor`: Noise spreading factor.
      - `poisson.NREM.poisson_weight`: Weight of Poisson input.
  - `weights`: Synaptic weight settings.
    - `weights.awake`: Settings for awake state.
      - `weights.awake.recurrent_weight`: Weight of recurrent synapses.
      - `weights.awake.inh_to_exc_weight`: Weight of inhibitory to excitatory synapses.
      - `weights.awake.exc_to_inh_weight`: Weight of excitatory to inhibitory synapses.
      - `weights.awake.inh_to_inh_weight`: Weight of inhibitory to inhibitory synapses.
      - `weights.awake.inter_pop_weight`: Weight of inter-population synapses.
    - `weights.NREM`: Settings for NREM state.
      - `weights.NREM.recurrent_weight`: Weight of recurrent synapses.
      - `weights.NREM.inh_to_exc_weight`: Weight of inhibitory to excitatory synapses.
      - `weights.NREM.exc_to_inh_weight`: Weight of excitatory to inhibitory synapses.
      - `weights.NREM.inh_to_inh_weight`: Weight of inhibitory to inhibitory synapses.
      - `weights.NREM.inter_pop_weight`: Weight of inter-population synapses.
  - `dc_exc`: Direct current settings for excitatory neurons.
    - `dc_exc.amplitudes`: Amplitudes for DC injection.
    - `dc_exc.delay`: Delay for DC injection.
    - `dc_exc.weight`: Weight of DC injection.
  - `dc_inh`: Direct current settings for inhibitory neurons.
    - `dc_inh.amplitude`: Amplitude for DC injection.
    - `dc_inh.delay`: Delay for DC injection.
    - `dc_inh.weight`: Weight of DC injection.

- **Usage**: Adjust these parameters to change the overall behavior of the simulation, such as the duration of the run or the specific neuron models used.

### 2. `basic_crop_and_plot.yaml`
- **Purpose**: Defines parameters for cropping data and generating plots.
- **Key Parameters**:
  - `crop`: Settings for cropping the data.
    - `crop.default_cropping`: Default cropping settings.
    - `crop.delay_ms`: Delay time in milliseconds for cropping.
    - `crop.duration_ms`: Duration in milliseconds for the data to be cropped.
  - `sampling`: Sampling settings for data analysis.
    - `sampling.spikes_sampling_window_ms`: Time window in milliseconds for spikes sampling.
    - `sampling.low_freq_sampling_window_ms`: Time window in milliseconds for low-frequency data sampling.
  - `plot`: Settings for generating plots.
    - `plot.save_plots_in_output_dir`: Boolean flag indicating whether plots should be saved in the output directory.
    - `plot.intermediate_plt_show`: Boolean flag to display intermediate plots.
    - `plot.excitatory_rastergram`: Settings for plotting excitatory rastergrams.
    - `plot.inhibitory_rastergram`: Settings for plotting inhibitory rastergrams.
    - `plot.instantaneous_firing_rate`: Settings for plotting instantaneous firing rates.
    - `plot.smoothed_spikes_firing_rate`: Settings for plotting smoothed firing rates based on spikes.
    - `plot.spectrogram_window`: Settings for the spectrogram window.
      - `plot.spectrogram_window.default_window`: Default window for spectrogram analysis.
      - `plot.spectrogram_window.custom_ms`: Custom window length in milliseconds for spectrogram analysis.
    - `plot.spike_activation_kernel`: Kernel settings for spike activation.
      - `plot.spike_activation_kernel.duration_ms`: Duration in milliseconds for the activation kernel.
      - `plot.spike_activation_kernel.std_dev_ms`: Standard deviation in milliseconds for the kernel.
      - `plot.spike_activation_kernel.plot`: Boolean flag to indicate whether to plot the activation kernel.
    - `plot.high_freq_spectrum`: Settings for plotting high-frequency spectrum.
      - `plot.high_freq_spectrum.max_plot_freq_Hz`: Maximum frequency in Hz for plotting the spectrum.
      - `plot.high_freq_spectrum.smoothing_length`: Smoothing length for the spectrum.
      - `plot.high_freq_spectrum.plot`: Boolean flag to indicate whether to plot the high-frequency spectrum.
    - `plot.high_freq_spectrogram`: Settings for high-frequency spectrogram plots.
      - `plot.high_freq_spectrogram.plot`: Boolean flag to indicate whether to plot the high-frequency spectrogram.
      - `plot.high_freq_spectrogram.max_plot_freq_Hz`: Maximum frequency in Hz for plotting the spectrogram.
    - `plot.tissue_response_kernel`: Kernel settings for tissue response.
      - `plot.tissue_response_kernel.duration_ms`: Duration in milliseconds for the tissue response kernel.
      - `plot.tissue_response_kernel.std_dev_ms`: Standard deviation in milliseconds for the kernel.
      - `plot.tissue_response_kernel.plot`: Boolean flag to indicate whether to plot the tissue response kernel.
    - `plot.tissue_response_rate`: Settings for plotting tissue response rates.
    - `plot.tissue_response_spectrum`: Settings for plotting tissue response spectrum.
      - `plot.tissue_response_spectrum.max_plot_freq_Hz`: Maximum frequency in Hz for the tissue response spectrum.
      - `plot.tissue_response_spectrum.smoothing_length`: Smoothing length for the spectrum.
      - `plot.tissue_response_spectrum.plot`: Boolean flag to indicate whether to plot the tissue response spectrum.
    - `plot.tissue_response_spectrogram`: Settings for plotting tissue response spectrogram.
      - `plot.tissue_response_spectrogram.plot`: Boolean flag to indicate whether to plot the tissue response spectrogram.
      - `plot.tissue_response_spectrogram.max_plot_freq_Hz`: Maximum frequency in Hz for plotting the spectrogram.

- **Usage**: Modify these settings to adjust which parts of the recorded data are visualized and how the plots are generated.

### 3. `basic_sim_and_recording_times.yaml`
- **Purpose**: Sets the timing parameters for simulations and recordings.
- **Key Parameters**:
  - `start_time`: The time at which the simulation starts.
  - `end_time`: The time at which the simulation ends.
  - `recording_intervals`: Specifies the intervals at which metrics such as firing rates are recorded.
  - `sampling_frequency`: Defines the frequency at which data should be sampled for analysis.
  - `recording_channels`: List of channels or populations to be recorded during the simulation.
- **Usage**: Use this file to define start and end times for simulations, as well as when to record specific metrics, such as firing rates.

### 4. `standard_AdEx_exc_neural_params.yaml` and `standard_AdEx_inh_neural_params.yaml`
- **Purpose**: These files define the parameters for excitatory and inhibitory neurons using the Adaptive Exponential (AdEx) model.
- **Key Parameters**:
  - `v_rest`: Resting membrane potential of the neuron.
  - `tau_m`: Membrane time constant.
  - `v_threshold`: Spike threshold potential.
  - `adaptation`: Parameters for spike frequency adaptation.
  - `refractory_period`: Duration of the refractory period after a spike.
  - `spike_amplitude`: Amplitude of the spike.
  - `reset_voltage`: Voltage to which the membrane potential is reset after a spike.
- **Usage**: Adjust these files to modify the behavior of excitatory and inhibitory neurons, such as their threshold potentials, time constants, and adaptation properties.

### 5. `Ca-AdEx_neural_params.yaml` and `MC-g0_like-exc_neural_params.yaml`
- **Purpose**: These files provide parameters for specific types of neurons, such as Calcium-modified AdEx or other custom neural models.
- **Key Parameters**:
  - `calcium_conductance`: Sets the conductance value for calcium ions, affecting neuron excitability.
  - `g0_exc`: Defines the baseline excitatory conductance.
  - `tau_calcium`: Time constant for calcium dynamics.
  - `calcium_threshold`: Threshold for calcium concentration to trigger certain behaviors.
  - `calcium_decay`: Defines the decay rate of calcium concentration.
- **Usage**: These files can be modified to experiment with different neural models and their properties within the network.

### 6. `basic_directories_and_list_of_yamls.yaml`
- **Purpose**: Specifies directories for saving outputs and lists the various YAML configuration files used in the simulation.
- **Key Parameters**:
  - `output_directory`: Directory path where results will be saved.
  - `yaml_files_list`: List of YAML files to be used during the simulation.
  - `overwrite_output`: Boolean value indicating whether to overwrite existing output files.
- **Usage**: Modify this file to set custom directories for storing results and to define which YAML files should be used for a specific run.

### 7. `tune_general_config.yaml`, `tune_sim_and_recording_times.yaml`, and `tune_crop_and_plot.yaml`
- **Purpose**: These files provide additional tuning options for the general configuration, simulation times, and plotting parameters.
- **Key Parameters**:
  - **`tune_general_config.yaml`**:
    - `learning_rate`: Sets the rate of synaptic plasticity during training phases.
    - `synaptic_scaling`: Enables or disables synaptic scaling to maintain homeostasis.
    - `input_scaling_factor`: Scales the strength of external inputs.
    - `background_noise`: Parameters for adding background noise to the network.
  - **`tune_sim_and_recording_times.yaml`**:
    - `adjust_start_time`: Fine-tune the start time for specific events within the simulation.
    - `adjust_end_time`: Fine-tune the end time for specific events.
    - `additional_recording_times`: Specifies additional time points for recording specific metrics.
  - **`tune_crop_and_plot.yaml`**:
    - `highlighted_regions`: Specifies regions of interest to highlight in the plots.
    - `plot_resolution`: Sets the resolution for plot outputs.
    - `color_scheme`: Defines the color scheme for visualizing different plots.
- **Usage**: Use these files to override default settings in the `basic_` configuration files, allowing for fine-tuning of specific simulation aspects without altering the base configurations.

## Troubleshooting

### Common Issues

- **Issue 1**: Jupyter Notebook fails to start

  - **Solution**: Ensure all dependencies are installed, and use Python 3.8 or later.

- **Issue 2**: NEST Simulator errors

  - **Solution**: Verify that NEST is installed correctly and that the environment variables are set.

## FAQ

**Q: How do I update the software?**A: Run the following command:

```bash
git pull origin main
```



