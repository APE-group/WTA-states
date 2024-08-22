# WTA-states
Exercises combining brain states and soft winner takes all

LOG OF MAIN EVENTS (from recent to past)

Release 2027-0822 solving issue #19 switch on plasticity and #23: Added capability of disconnecting intra assemblies connections at programmable time, and what target exc pop to be disconnected

Release 2027-0820b config and plots can be saved in directories that can be specified using the file "directories_and_list_of_yamls.yaml".
Use there the "overwrite_dir" file name if you want to send last runs always to the same directory.

Release 2027-0820 Solution of issue #16. It is now possible to specify a list of events during the simulation. The first supported action is the saving of intra assembly synaptic matrices at times that can be defined in the basic_ and tune_sim_and_recording_times.yaml files

Release 2024-0817 It Solves both issue #11 and issue #12. Now the main launches a single activity analysis function. 
All parameters for the analysis and the plot of the activity can be set in the basic_ and tune_ crop_and_plot.yaml files.

Commit 2024-0817 Solving issue #12 create crop_and_plot yaml files

Release 2024-0816 Solving issue #9: "encapsulate the first part of the main in a nest_parameters_preparation() function"

Release 2024-0808 Solving issue #5: "Add basic contextual signal reaching a single assembly with a start and stop time"it is the product of branch: "#5-add-basic-contextual-signal..."

Release 2024-0807 introduced multicompartment support solving issue#3: "Support for runs with both single and multi-compartment neurons"

Release 2024-0510 Added capability of spectral analysis with error bands, by pull request #1 that solves issue #2: "Need for spectral analysis with error band and detection of maxima"

2024-0430: novel features:
- all parameters stored in reference yaml files
- all parameters can be tuned reading additional "tune_..." yaml files 
- nest network creation and simulation in separate file

2024-0428 Initial commit 
