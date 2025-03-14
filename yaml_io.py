import yaml
import os
import shutil

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

def read_crop_and_plot_yaml(directories_and_list_of_yamls, verbose):
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
        
    #add info about possible absolute path to store plots
    current_dir_name = os.getcwd() 
    relative_output_dir_name =\
      directories_and_list_of_yamls['relative_output_dir']
    absolute_plot_dir_path = current_dir_name + "/" + relative_output_dir_name  
    plot_params['absolute_plot_dir_path'] = absolute_plot_dir_path

    return crop_params, plot_params, sampling_params

def read_basic_directories_and_list_of_yamls(verbose):
    print("in read_directories_yaml: verbose mode is", verbose)
    directories_and_yaml_files_name ='basic_directories_and_list_of_yamls.yaml'
    with open(directories_and_yaml_files_name, 'r') as file:
        directories_and_list_of_yamls = yaml.safe_load(file)
    return directories_and_list_of_yamls

def clean_directory(dir_path):
    """
    Removes all files and subdirectories inside the given directory.

    dir_path: The path to the directory you want to clean.
    """
    # List everything in the directory
    for item in os.listdir(dir_path):
        item_path = os.path.join(dir_path, item)
        
        # Check whether it's a file, symlink, or directory
        if os.path.isfile(item_path) or os.path.islink(item_path):
            os.remove(item_path)   # remove file or symlink
        elif os.path.isdir(item_path):
            shutil.rmtree(item_path)  # remove entire directory tree

# Example usage:
# clean_directory("/path/to/your/directory")

def copy_yamls_in_output_dir(directories_and_list_of_yamls,verbose):
    # Define the directory and file paths
    current_dir_name = os.getcwd() 
    relative_output_dir_name=\
      directories_and_list_of_yamls['relative_output_dir']
    
    if verbose:
        print("in copy_yamls_in_output_directory: current execution directory:", current_dir_name)
        print("in copy_yamls_in_output_directory: relative_output_dir:", relative_output_dir_name)
    output_dir_name = current_dir_name + "/" + relative_output_dir_name
    if verbose:
        print("in copy_yamls_in_output_directory: output dir:", output_dir_name)    
  
    # Check if the directory exists
    if not os.path.exists(output_dir_name):
        os.makedirs(output_dir_name)
        print(f"Directory '{output_dir_name}' created.")
    else:
        if relative_output_dir_name != 'overwrite_dir':
            response = input(f"Directory '{output_dir_name}' already exists. Do you want to continue? (y/n): ")
            if response.lower() != 'y':
                print("Operation aborted.")
                exit()
    clean_directory(output_dir_name)
    file_names = directories_and_list_of_yamls['config_yaml_files']
    # Copy the files to the directory
    
    for file_name in file_names: 
        input_file_name = current_dir_name + "/" + file_name
        destination_file_name = output_dir_name + "/" + file_name
        shutil.copy(input_file_name, destination_file_name)
        if verbose:
            print(f"File '{input_file_name}' copied to '{destination_file_name}'.")

    return output_dir_name


def copy_neu_params_yamls_in_output_dir(output_dir_name, config, verbose):
    assert(os.path.exists(output_dir_name))
    current_dir_name = os.getcwd()
    
    #copying the yaml describing inhibitory neurons
    input_file_name = current_dir_name + "/" + config['inh_neu_params_filename']
    destination_file_name = output_dir_name + "/" + config['inh_neu_params_filename']
    shutil.copy(input_file_name, destination_file_name)
    if verbose:
        print(f"File '{input_file_name}' copied to '{destination_file_name}'.")

    #copying the yaml describing excitatory neurons
    input_file_name = current_dir_name + "/" + config['exc_neu_params_filename']
    destination_file_name = output_dir_name + "/" + config['exc_neu_params_filename']
    shutil.copy(input_file_name, destination_file_name)
    if verbose:
        print(f"File '{input_file_name}' copied to '{destination_file_name}'.")
    return
            
        
        
        