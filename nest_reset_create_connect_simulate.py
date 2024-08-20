import nest
from cm_neuron import create_cm_neuron
from synaptic_io import *


def nest_reset_create_connect_simulate(nest_pms, num_threads, verbose):
    sim_completed=False
    sim_pms=nest_pms["sim_pms"]
    nest.ResetKernel()

    nest.SetKernelStatus({"resolution": sim_pms["resolution_ms"]})
    nest.SetKernelStatus({'local_num_threads':num_threads})
    
    num_exc_neu_per_pop = nest_pms["network"]["num_exc_neu_per_pop"]
    num_exc_pop = nest_pms["network"]["num_exc_pop"]

    use_single_compartment_environment=nest_pms["use_single_compartment_environment"]
    print("IN nest_reset_create_connect_simulate: use_single_compartment_environment =", use_single_compartment_environment)
    exc_neu_params=nest_pms['exc_neu_params'] 

    if use_single_compartment_environment:
        # Create excitatory neuron populations
        neurons = [
            nest.Create(
                "aeif_cond_alpha",
                num_exc_neu_per_pop,
                params=exc_neu_params['equation_params'],
            )
            for _ in range(num_exc_pop)
        ]
    else:
        # Create excitatory Ca-AdEx neuron populations
        nest.CopyModel('static_synapse', 'REFRACT_syn')
        nest.CopyModel('static_synapse', 'ADAPT_syn')
        nest.CopyModel('static_synapse', 'BAP_syn')
        
        neurons = [
            create_cm_neuron(num_exc_neu_per_pop, params=exc_neu_params['equation_params']) for _ in range(num_exc_pop)
        ]

    # Create inhibitory neurons
    inh_neu_params = nest_pms["inh_neu_params"]
    num_inh_neu = nest_pms["network"]["num_inh_neu"]
    inh_neurons = nest.Create(
        "aeif_cond_alpha",
        num_inh_neu,
        params=inh_neu_params,
    )

    # Spike recorders for excitatory and inhibitory populations
    recording_pms=nest_pms["recording_pms"]
    spike_recorders = [nest.Create("spike_recorder",
        {"start": recording_pms["start_ms"], "stop": recording_pms["stop_ms"]}) 
                   for _ in range(num_exc_pop)]
    inh_spike_recorder = nest.Create("spike_recorder", 
        {"start": recording_pms["start_ms"], "stop": recording_pms["stop_ms"]})

    # Connect spike detectors
    for i in range(num_exc_pop):
        nest.Connect(neurons[i], spike_recorders[i])
    nest.Connect(inh_neurons, inh_spike_recorder)

    # Connection specifications
    exc_t_ref_ms=nest_pms["exc_t_ref_ms"]
    min_syn_delay_ms=nest_pms["network"]["min_syn_delay_ms"]
    use_recurrency=nest_pms["network"]["use_exc_recurrency"]
    recurrent_weight = nest_pms["network"]["weights"]["exc_to_exc"] 
    inh_to_exc_weight = nest_pms["network"]["weights"]["inh_to_exc_weight"]
    exc_to_inh_weight = nest_pms["network"]["weights"]["exc_to_inh_weight"] 
    inh_to_inh_weight = nest_pms["network"]["weights"]["inh_to_inh_weight"]
    
    if use_recurrency:
        exc_to_exc_delay_ms = min_syn_delay_ms + exc_t_ref_ms + 1.0
        inh_to_inh_delay_ms = min_syn_delay_ms + exc_t_ref_ms + 1.0
        conn_spec_dict = {"rule": "all_to_all", "allow_autapses": False}
        if(use_single_compartment_environment):
            for i in range(num_exc_pop):
                nest.Connect(neurons[i], neurons[i], conn_spec_dict,\
                    syn_spec={"weight": recurrent_weight, "delay": exc_to_exc_delay_ms})
        else:
            for i in range(num_exc_pop):
                nest.Connect(neurons[i], neurons[i], conn_spec_dict,\
                    syn_spec={"weight": recurrent_weight, "delay": exc_to_exc_delay_ms, 'receptor_type': 9})
        nest.Connect(inh_neurons, inh_neurons, conn_spec_dict,\
                    syn_spec={"weight": inh_to_inh_weight, "delay": inh_to_inh_delay_ms})
    
    # Connect inhibitory neurons to all excitatory neurons and vice versa
    inh_to_exc_delay_ms = min_syn_delay_ms + exc_t_ref_ms + 0.55
    exc_to_inh_delay_ms = min_syn_delay_ms + 0.5
    for pop in neurons:
        if use_single_compartment_environment:        
            nest.Connect(inh_neurons, pop, {"rule": "all_to_all"},\
                syn_spec={"weight": inh_to_exc_weight, "delay": inh_to_exc_delay_ms})
            nest.Connect(pop, inh_neurons, {"rule": "all_to_all"},\
                syn_spec={"weight": exc_to_inh_weight, "delay": exc_to_inh_delay_ms})
        else:
            nest.Connect(inh_neurons, pop, {"rule": "all_to_all"},\
                syn_spec={"weight": -inh_to_exc_weight, "delay": inh_to_exc_delay_ms, 'receptor_type': 10})
            nest.Connect(pop, inh_neurons, {"rule": "all_to_all"},\
                syn_spec={"weight": exc_to_inh_weight, "delay": exc_to_inh_delay_ms})

    # Create and connect Poisson generators if enabled
    use_poisson_generators=nest_pms["network"]["use_poisson_generators"]
    if use_poisson_generators:
        num_poisson_generators=nest_pms["network"]["num_poisson_generators"]
        poisson_rate=nest_pms["poisson"]["poisson_rate"]
        pgs = nest.Create("poisson_generator", num_poisson_generators, params={"rate": poisson_rate})
        poisson_weight=nest_pms["poisson"]["poisson_weight"]
        for i in range(num_exc_pop):
            if use_single_compartment_environment:
                nest.Connect(pgs, neurons[i],\
                             syn_spec={"weight": poisson_weight, "delay": 1.0})
            else:
                nest.Connect(pgs, neurons[i], syn_spec={"weight": poisson_weight, "delay": 1.0, 'receptor_type': 9})

    # Add contextual poisson signal if configured
    if 'contextual_poisson' in nest_pms and \
            nest_pms['brain_state'] in nest_pms['contextual_poisson']: 
        # Extract the upper-level configuration parameters
        contextual_poisson_config = nest_pms['contextual_poisson'][nest_pms['brain_state']]
        
        # Ensure all required keys are present in the configuration
        if all(key in contextual_poisson_config for key in ['events', 'spreading_factor', 'basic_rate', 'poisson_weight']):
            events = contextual_poisson_config['events']
            spreading_factor = contextual_poisson_config['spreading_factor']
            basic_rate = contextual_poisson_config['basic_rate']
            poisson_weight = contextual_poisson_config['poisson_weight']
    
            # Iterate through each event in the events list
            for event in events:
                target_pop = event['target_population']
                start_time = event['start_time_ms']
                stop_time = event['stop_time_ms']
                
                # Create a Poisson generator for the event
                contextual_poisson_gen = nest.Create('poisson_generator', 1, {'rate': basic_rate})
                
                # Set the start and stop times for the generator
                nest.SetStatus(contextual_poisson_gen, {'start': start_time, 'stop': stop_time})
    
                # Define the synapse specifications
                syn_spec = {'weight': poisson_weight * spreading_factor, "delay": 1.0}
    
                # Connect the generator to the target population
                if use_single_compartment_environment:             
                    nest.Connect(contextual_poisson_gen, neurons[target_pop], syn_spec=syn_spec)
                else:
                    syn_spec.update({'receptor_type': 9})
                    nest.Connect(contextual_poisson_gen, neurons[target_pop], syn_spec=syn_spec)


    # DC current injection for all neurons if enabled
    use_dc_exc_injectors=nest_pms["use_dc_exc_injectors"]
    if use_dc_exc_injectors:
        dc_exc_start_ms=nest_pms["dc_exc"]["start_ms"]
        dc_exc_stop_ms=nest_pms["dc_exc"]["stop_ms"]
        dc_exc_amplitudes=nest_pms["dc_exc"]["dc_exc_amplitudes"]
        dc_exc_weight=nest_pms["dc_exc"]["dc_exc_weight"]
        dc_exc_delay_ms=nest_pms["dc_exc"]["dc_exc_delay_ms"]
        dcgs = [nest.Create("dc_generator", 1,\
                {"start": dc_exc_start_ms, "stop": dc_exc_stop_ms, "amplitude": amp})\
                for amp in dc_exc_amplitudes]
        for i in range(num_exc_pop):
            if use_single_compartment_environment:
                nest.Connect(dcgs[i], neurons[i],\
                    syn_spec={'weight': dc_exc_weight, 'delay': dc_exc_delay_ms})
            else:
                #assert("MC neuron not yet implemented")
                nest.Connect(dcgs[i], neurons[i], syn_spec={'weight': dc_exc_weight, 'delay': dc_exc_delay, 'receptor_type': 0})

    # Specific DC current injection for inhibitory neurons
    use_dc_inh_injector=nest_pms["use_dc_inh_injector"]
    if use_dc_inh_injector:
        dc_inh_start_ms = nest_pms["dc_inh"]["start_ms"]
        dc_inh_stop_ms = nest_pms["dc_inh"]["stop_ms"]
        dc_inh_amplitude = nest_pms["dc_inh"]["dc_inh_amplitude"]
        dc_inh_weight = nest_pms["dc_inh"]["weight"]
        dc_inh_delay_ms = nest_pms["dc_inh"]["delay_ms"]
        
        dc_inh_generator = nest.Create("dc_generator", 1,\
                {"start": dc_inh_start_ms, "stop": dc_inh_stop_ms, 
                 "amplitude": dc_inh_amplitude})
        nest.Connect(dc_inh_generator, inh_neurons,\
                     syn_spec={'weight': dc_inh_weight, 'delay': dc_inh_delay_ms})

    def perform_event_action(requested_action,verbose):
        if verbose:
            print("in perform_event_action: requested action:", requested_action) 
        if requested_action['kind']=='store_intra_assembly_syn':
            store_intra_assembly_syn(requested_action['kind'],requested_action['file_name'],verbose)
        else:
            print("in perform_event_action: requested action:",requested_action['kind'],"NOT SUPPORTED")
            assert False
        return

    def store_intra_assembly_syn(requested_action, file_name, verbose):
        assert(requested_action=='store_intra_assembly_syn')
        if(verbose):
            print("in store_intra_assembly_synapses:")
            print("will save to file:", file_name)
        array_of_dicts = []
        for i in range(num_exc_pop):
            array_of_dicts.append(\
                {"exc_pop_index": i,\
                "intra_pop_connections": nest.GetConnections(source=neurons[i],target=neurons[i]).get(
                ("source", "target", "weight"), output="pandas")})
        with open(file_name, "wb") as f:
            pickle.dump(array_of_dicts, f, pickle.HIGHEST_PROTOCOL)
        
        #print("in store_intra_assembly_synapses: ACTION NOT YET IMPLEMENTED")
        return
    
    # Run the simulation
    verbose=True
    if verbose:
        print("in nest_..._simulate: just before nest.Simulate")
        print("nest_pms['events_pms']", nest_pms['events_pms'])
    tot_simulated_time_ms=0.0
    for item in nest_pms['events_pms']:
        if verbose:
            print("ACTION",item['action'])
        #performing the action
        perform_event_action(item['action'],verbose)   

        if verbose:       
            print("expected NEXT SIM PERIOD ms",item['sim_period_ms']) 
        next_sim_period_ms = item['sim_period_ms']
        if next_sim_period_ms + tot_simulated_time_ms > sim_pms["stop_ms"]:
            next_sim_period_ms = sim_pms["stop_ms"] - tot_simulated_time_ms
            if verbose:
                print("reduced last sim time ms to", next_sim_period_ms) 
        if verbose:
            print("launching next period of ", next_sim_period_ms, "ms of simulation")
        #launching simulation period
        nest.Simulate(next_sim_period_ms)
        tot_simulated_time_ms += next_sim_period_ms
        if verbose:
            print("performed until now:",tot_simulated_time_ms, "TOTAL ms of simulation")
    if tot_simulated_time_ms < sim_pms["stop_ms"]:
        next_sim_period_ms = sim_pms["stop_ms"] - tot_simulated_time_ms
        if verbose:
            print("ADDING ", next_sim_period_ms, "ms sim to complete expected" ,sim_pms["stop_ms"], "TOTAL ms of simulation")
        #possible last simulation period
        nest.Simulate(next_sim_period_ms)
        tot_simulated_time_ms += next_sim_period_ms
        if verbose:
            print("performed until now:",tot_simulated_time_ms, "TOTAL ms of simulation")
        
    sim_completed=True
    return sim_completed, spike_recorders, inh_spike_recorder
