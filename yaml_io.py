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
    with open("sim_and_recording_times.yaml", 'r') as file:
        times = yaml.safe_load(file)
    if verbose:
        print("ORIGINAL YAML", times)
        print("----")
    with open("tune_sim_and_recording_times.yaml", 'r') as file:
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
    # Load standard network parameters from YAML file
    with open("general_config.yaml", 'r') as file:
        config = yaml.safe_load(file)
    if verbose:
        print("ORIGINAL CONFIG", config)
        print("----")     
    with open("tune_general_config.yaml", 'r') as file:
        tune_config = yaml.safe_load(file)
    
    # Use the recursive update function
    recursive_update(config, tune_config)

    if verbose:
        print("TUNED CONFIG", config)
        print("----")     
    return config

def read_neural_parameters(is_standard_nest,is_verbose):
    # Load neural parameters from a YAML file
    if is_standard_nest:
        neural_param_file = "standard_nest_neural_params.yaml" 
    else:
        print("MC not yet implemented")
        assert False
    with open(neural_param_file, 'r') as file:
        params = yaml.safe_load(file)
    exc_pms = params['excitatory']
    inh_pms = params['inhibitory']
    return exc_pms,inh_pms