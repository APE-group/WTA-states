import pickle

import numpy as np

def load_intra_assembly_syn(file_name, verbose):
    if verbose:
        print("in load_intra_assembly_syns:")
        print("reading from file:", file_name)
    
    # Load the stored synaptic matrices
    with open(file_name, "rb") as f:
        array_of_dicts = pickle.load(f)
    
    # Calculate and print the mean and standard deviation for each population
    for pop_dict in array_of_dicts:
        exc_pop_index = pop_dict['exc_pop_index']
        synaptic_weights = pop_dict['intra_pop_connections']['weight'].values  #         Extract the weights as a numpy array
        mean_weight = np.mean(synaptic_weights)
        std_weight = np.std(synaptic_weights)
        
        print(f"Population {exc_pop_index}:")
        print(f"  Mean Weight: {mean_weight:.4f}")
        print(f"  Standard Deviation: {std_weight:.4f}")
    
    return array_of_dicts
