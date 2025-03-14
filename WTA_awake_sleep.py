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
from nest_reset_create_connect_simulate import *
from activity_analysis import *
from synaptic_io import *
from synaptic_analysis import *

def WTA_awake_sleep():

    is_verbose = False
    #copy configuration yamls in specified output directory
    directories_and_list_of_yamls =\
        read_basic_directories_and_list_of_yamls(is_verbose)
    copy_yamls_in_output_dir(directories_and_list_of_yamls, is_verbose)

    is_verbose = True
    #total sim, resolution and recording times
    times = read_sim_and_recording_times_yaml(is_verbose)

    #read general network parameters
    config = read_general_config_yaml(is_verbose)

    # copy neural params files in specified output directory
    output_dir_name = directories_and_list_of_yamls['relative_output_dir']
    copy_neu_params_yamls_in_output_dir(output_dir_name, config, is_verbose)

    #prepare all simulation parameters
    nest_pms={}
    nest_pms = nest_parameters_preparation(times, config, is_verbose, nest_pms)
    print("nest_pms",nest_pms)

    NEST_version = nest.__version__
    if NEST_version == "3.7.0" and nest_pms["use_single_compartment_environment"]==False:
        print("ASSERTION ERROR: Ca-AdEx multi-compartment neuron not supported by this NEST version", NEST_version)
        assert(False)

    is_verbose=False
    num_threads=4
    sim_completed, spike_recorders, inh_spike_recorder, multimeters, list_of_syn_matrix_file_names = nest_reset_create_connect_simulate(nest_pms,num_threads, is_verbose)
    print("sim_completed", sim_completed)

    d_inh = nest.GetStatus(inh_spike_recorder, "events")[0]

    #before analysis, preliminary sim look 
    is_verbose=True
    preliminary_sim_look(is_verbose, nest_pms, spike_recorders, inh_spike_recorder, nest_pms["recording_pms"])   

    is_verbose=False
    #here we prepare all the parameters for the following analysis and print
    crop_pms, plot_pms, sampling_pms, analysis_pms = prepare_crop_plot_sampling_activityAnalysis_parameters(directories_and_list_of_yamls, nest_pms, is_verbose)

    # Assuming spike_recorders is a list of spike recorder IDs previously created in your NEST simulation
    cropped_events = crop_events_from_spike_recorders(crop_pms, spike_recorders)
    cropped_inh_events = crop_inh_events(crop_pms, inh_spike_recorder)

    is_verbose=False
    #launches all analysis
    #produces from both spike-like waveforms and tissue-like responses
    #PLEASE use the basic_ and tune_crop_and_plot.yaml to select plots and parameters 
    produce_rastegrams_rates_spectra_spectrograms(directories_and_list_of_yamls,nest_pms, crop_pms, plot_pms, analysis_pms, cropped_events, cropped_inh_events, is_verbose)

    verbose=False
    # loading and printing of info from synaptic matrices
    for file_name in list_of_syn_matrix_file_names:
        array_of_dicts=load_syn(file_name, verbose) 

        syn_analysis(list_of_syn_matrix_file_names)

    return



