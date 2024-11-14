import numpy as np
import matplotlib.pyplot as plt
import os
import pickle


def syn_analysis(list_of_syn_matrix_file_names):
    """
    Analyzes and plots the temporal evolution of synaptic matrices.

    Parameters:
        list_of_syn_matrix_file_names (list): A list of file names containing synaptic matrices saved at different times.
    """
    time_points = []  # To store the times (in ms) when the synaptic matrices were saved
    synaptic_data = {}  # To store synaptic weights for each source-target pair over time

    # Iterate through each file to extract data
    for file_name in list_of_syn_matrix_file_names:
        if not os.path.exists(file_name):
            print(f"File {file_name} does not exist. Skipping.")
            continue

        # Extracting the synaptic information from the pickle file
        try:
            with open(file_name, "rb") as f:
                array_of_dicts = pickle.load(f)
        except Exception as e:
            print(f"Error reading pickle file {file_name}: {e}. Skipping.")
            continue

        # Iterate through the list of synaptic data dictionaries
        for pop_dict in array_of_dicts:
            time_ms = pop_dict.get('time (ms)', None)
            if time_ms is None:
                print(f"Time information not found in file {file_name}. Skipping entry.")
                continue

            source_target = pop_dict.get('conn_index', None)  # Assuming 'conn_index' identifies the source-target pair
            synaptic_weights = pop_dict['connections']['weight'].values  # Extract the weights as a numpy array

            # Add the time point to the list if it's not already added
            if time_ms not in time_points:
                time_points.append(time_ms)

            # Store synaptic weights for each source-target pair
            if source_target not in synaptic_data:
                synaptic_data[source_target] = {}
            if time_ms not in synaptic_data[source_target]:
                synaptic_data[source_target][time_ms] = []
            synaptic_data[source_target][time_ms].extend(synaptic_weights)

    # Sort time points and corresponding data for proper plotting
    time_points = sorted(time_points)

    # Plotting the temporal evolution of synaptic weights
    plt.figure(figsize=(10, 6))

    markers = ['o', 's', '^', 'D', 'v', '>', '<', 'p', '*']  # Different markers for different source-target pairs
    marker_index = 0

    for source_target, weight_dict in synaptic_data.items():
        mean_weights = []
        std_weights = []

        # Calculate mean and standard deviation of synaptic weights across time
        for time in time_points:
            weights_at_time = weight_dict.get(time, [])
            if weights_at_time:
                mean_weights.append(np.mean(weights_at_time))
                std_weights.append(np.std(weights_at_time))
            else:
                mean_weights.append(np.nan)
                std_weights.append(np.nan)

        # Plot with error bars and different markers for each source-target pair
        marker = markers[marker_index % len(markers)]
        plt.errorbar(time_points, mean_weights, yerr=std_weights, label=f'{source_target} (marker: {marker})', capsize=5, marker=marker, linestyle='-')
        marker_index += 1

    plt.xlabel('Time (ms)')
    plt.ylabel('Weight')
    plt.title('Temporal Evolution of Synaptic Weights')
    plt.legend()
    plt.grid(True)
    plt.tight_layout()
    plt.show()

# Example usage
# syn_analysis(["syn_matrix_100ms.pkl", "syn_matrix_200ms.pkl", "syn_matrix_300ms.pkl"])

"""


import numpy as np
import matplotlib.pyplot as plt
import os
import pickle


def syn_analysis(list_of_syn_matrix_file_names):
    
"""
"""    
    Analyzes and plots the temporal evolution of synaptic matrices.

    Parameters:
        list_of_syn_matrix_file_names (list): A list of file names containing synaptic matrices saved at different times.
"""
"""
    time_points = []  # To store the times (in ms) when the synaptic matrices were saved
    synaptic_data = {}  # To store synaptic weights for each assembly over time

    # Iterate through each file to extract data
    for file_name in list_of_syn_matrix_file_names:
        if not os.path.exists(file_name):
            print(f"File {file_name} does not exist. Skipping.")
            continue

        # Extracting the synaptic information from the pickle file
        try:
            with open(file_name, "rb") as f:
                array_of_dicts = pickle.load(f)
        except Exception as e:
            print(f"Error reading pickle file {file_name}: {e}. Skipping.")
            continue

        # Iterate through the list of synaptic data dictionaries
        for pop_dict in array_of_dicts:
            time_ms = pop_dict.get('time (ms)', None)
            if time_ms is None:
                print(f"Time information not found in file {file_name}. Skipping entry.")
                continue

            synapse_model = pop_dict.get('synapse_model', None)
            synaptic_weights = pop_dict['connections']['weight'].values  # Extract the weights as a numpy array

            # Add the time point to the list if it's not already added
            if time_ms not in time_points:
                time_points.append(time_ms)

            # Store synaptic weights for each synapse model
            if synapse_model not in synaptic_data:
                synaptic_data[synapse_model] = {}
            if time_ms not in synaptic_data[synapse_model]:
                synaptic_data[synapse_model][time_ms] = []
            synaptic_data[synapse_model][time_ms].extend(synaptic_weights)

    # Sort time points and corresponding data for proper plotting
    time_points = sorted(time_points)

    # Plotting the temporal evolution of synaptic weights
    plt.figure(figsize=(10, 6))

    for synapse_model, weight_dict in synaptic_data.items():
        mean_weights = []
        std_weights = []

        # Calculate mean and standard deviation of synaptic weights across time
        for time in time_points:
            weights_at_time = weight_dict.get(time, [])
            if weights_at_time:
                mean_weights.append(np.mean(weights_at_time))
                std_weights.append(np.std(weights_at_time))
            else:
                mean_weights.append(np.nan)
                std_weights.append(np.nan)

        # Plot with error bars
        plt.errorbar(time_points, mean_weights, yerr=std_weights, label=f'{synapse_model}', capsize=5)

    plt.xlabel('Time (ms)')
    plt.ylabel('Weight')
    plt.title('Temporal Evolution of Synaptic Weights')
    plt.legend()
    plt.grid(True)
    plt.tight_layout()
    plt.show()

# Example usage
# syn_analysis(["syn_matrix_100ms.pkl", "syn_matrix_200ms.pkl", "syn_matrix_300ms.pkl"])
"""


