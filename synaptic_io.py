import pickle

import numpy as np

def load_syn(file_name, verbose):
    if verbose:
        print("in load_syn:")
        print("reading from file:", file_name)
    
    # Load the stored synaptic matrices
    with open(file_name, "rb") as f:
        array_of_dicts = pickle.load(f)
    print("------------------")
    print("SYNAPTIC EVOLUTION")
    # Calculate and print the mean and standard deviation for each population
    for pop_dict in array_of_dicts:
        time_ms = pop_dict['time (ms)']
        index = pop_dict['conn_index']
        synapse_model = pop_dict['synapse_model']
        synaptic_weights = pop_dict['connections']['weight'].values  #         Extract the weights as a numpy array
        mean_weight = np.mean(synaptic_weights)
        std_weight = np.std(synaptic_weights)
        print("synapse_model",synapse_model)
        print(f"  Mean Weight: {mean_weight:.4f}")
        print(f"  Standard Deviation: {std_weight:.4f}")
    
    return array_of_dicts
