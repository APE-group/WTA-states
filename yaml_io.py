import yaml

def recursive_update(orig_dict, new_dict):
    """
    Recursively updates the original dictionary with values from the new dictionary.
    This function merges nested dictionaries instead of overwriting the entire value,
    preserving existing keys not present in the new dictionary.
    """
    for key, value in new_dict.items():
        if isinstance(value, dict) and key in orig_dict and isinstance(orig_dict[key], dict):
            recursive_update(orig_dict[key], value)
        else:
            orig_dict[key] = value

def read_sim_and_recording_times_yaml(verbose):
    print("in read_sim_and_recording_times_yaml: verbose mode is", verbose)
    default_sim_and_recording_times_file_name = "basic_sim_and_recording_times.yaml"
    tune_sim_and_recording_times_file_name = "tune_sim_and_recording_times.yaml"
    if verbose:
        print("IN read_sim_and_recording_times_yaml:") 
        print("opening file:", default_sim_and_recording_times_file_name)
    with open(default_sim_and_recording_times_file_name, 'r') as file:
        times = yaml.safe_load(file)
    if verbose:
        print("DEFAULT CONFIG (possibly to be TUNED):", times)
        print("----")
    if verbose:
        print("IN read_sim_and_recording_times_yaml:") 
        print("opening file:", tune_sim_and_recording_times_file_name)
    with open(tune_sim_and_recording_times_file_name, 'r') as file:
        tune_times = yaml.safe_load(file)
        
    # Use the recursive update function
    recursive_update(times, tune_times)  

    assert(times['sim_pms']['stop_ms'] >= 0.0 ) 
    assert(times['recording_pms']['start_ms'] <= times['sim_pms']['stop_ms'])
    assert(times['recording_pms']['start_ms'] <= times['recording_pms']['stop_ms']) 
    assert(times['recording_pms']['stop_ms'] <= times['sim_pms']['stop_ms'])  
    if verbose:
        print("TUNED CONFIG", times)
        print("----")  
    times["recording_pms"]["duration_ms"]=\
        times["recording_pms"]["stop_ms"]-\
        times["recording_pms"]["start_ms"]
    return times

def read_general_config_yaml(verbose):
    print("in read_general_config_yaml: verbose mode is", verbose)
    # Load standard network parameters from YAML file
    general_config_file_name = "basic_general_config.yaml"
    tune_config_file_name = "tune_general_config.yaml"
    if verbose:
        print("read_general_config_yaml:")
        print("read_general_config_yaml: opening file:",general_config_file_name)  
    with open(general_config_file_name, 'r') as file:
        config = yaml.safe_load(file)
    if verbose:
        print("read_general_config_yaml: starting config (later tuned by this routine)", config)
        print("----")
        print("read_general_config_yaml: opening file:",tune_config_file_name)  
    with open(tune_config_file_name, 'r') as file:
        tune_config = yaml.safe_load(file)
    
    # Use the recursive update function
    recursive_update(config, tune_config)

    if verbose:
        print("read_general_config_yaml: TUNED CONFIG", config)
    return config

def read_neural_parameters(neural_param_file_name, verbose):
    print("in read_neural_parameters: verbose mode is", verbose)
    # Load neural parameters from a YAML file
    if verbose:
        print("read_neural_parameters opening file:", neural_param_file_name) 
    with open(neural_param_file_name, 'r') as file:
        neural_params = yaml.safe_load(file)
    if verbose:    
        print("read_neural_parameters loaded the following parameters:", neural_params)
    return neural_params

def read_crop_and_plot_yaml(verbose):
    print("in read_crop_and_plot_yaml: verbose mode is", verbose)
    # Load standard crop and plot parameters from YAML file
    general_crop_and_plot_file_name = "basic_crop_and_plot.yaml"
    tune_crop_and_plot_file_name = "tune_crop_and_plot.yaml"
    if verbose:
        print("in read__crop_and_plot_yaml: opening file:",general_crop_and_plot_file_name)  
    with open(general_crop_and_plot_file_name, 'r') as file:
        crop_and_plot_params = yaml.safe_load(file)
    if verbose:
        print("in read__crop_and_plot_yaml: starting config (later tuned by this routine)", crop_and_plot_params)
        print("----")
        print("in read__crop_and_plot_yaml:: opening file:",tune_crop_and_plot_file_name)  
    with open(tune_crop_and_plot_file_name, 'r') as file:
        tune_crop_and_plot_params = yaml.safe_load(file)
    
    # Use the recursive update function
    recursive_update(crop_and_plot_params, tune_crop_and_plot_params)

    if verbose:
        print("in read__crop_and_plot_yaml: TUNED CONFIG", crop_and_plot_params)
        print("----") 
    crop_params=crop_and_plot_params['crop']
    if verbose:
        print("crop_params",crop_params)
    plot_params=crop_and_plot_params['plot']
    if verbose:
        print("plot_params",plot_params)
    sampling_params=crop_and_plot_params['sampling']
    if verbose:
        print("sampling_params",sampling_params) 
        
    return crop_params, plot_params, sampling_params