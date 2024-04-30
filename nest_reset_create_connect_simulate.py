import nest
def nest_reset_create_connect_simulate(nest_pms):
    sim_completed=False
    sim_pms=nest_pms["sim_pms"]
    #Creation of actual simulation starts herefrom
    nest.ResetKernel()

    nest.SetKernelStatus({"resolution": sim_pms["resolution_ms"]})
    
    num_exc_neu_per_pop = nest_pms["network"]["num_exc_neu_per_pop"]
    num_exc_pop = nest_pms["network"]["num_exc_pop"]

    use_standard_exc_nest_neuron=nest_pms["use_standard_exc_nest_neuron"]
    if use_standard_exc_nest_neuron:
        exc_neu_params=nest_pms["exc_neu_params"]
        # Create excitatory neuron populations
        neurons = [
            nest.Create(
                "aeif_cond_alpha",
                num_exc_neu_per_pop,
                params=exc_neu_params,
            )
            for _ in range(num_exc_pop)
        ]
    else:
        assert("MC exc neuron not yet supported")

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
        if(use_standard_exc_nest_neuron):
            for i in range(num_exc_pop):
                nest.Connect(neurons[i], neurons[i], conn_spec_dict,\
                    syn_spec={"weight": recurrent_weight, "delay": exc_to_exc_delay_ms})
        else:
            assert("MC neuron not yet implemented")
        nest.Connect(inh_neurons, inh_neurons, conn_spec_dict,\
                    syn_spec={"weight": inh_to_inh_weight, "delay": inh_to_inh_delay_ms})
    
    # Connect inhibitory neurons to all excitatory neurons and vice versa
    inh_to_exc_delay_ms = min_syn_delay_ms + exc_t_ref_ms + 0.55
    exc_to_inh_delay_ms = min_syn_delay_ms + 0.5
    for pop in neurons:
        if use_standard_exc_nest_neuron:        
            nest.Connect(inh_neurons, pop, {"rule": "all_to_all"},\
                syn_spec={"weight": inh_to_exc_weight, "delay": inh_to_exc_delay_ms})
            nest.Connect(pop, inh_neurons, {"rule": "all_to_all"},\
                syn_spec={"weight": exc_to_inh_weight, "delay": exc_to_inh_delay_ms})
        else:
            assert("MC neuron not yet implemented")

    # Create and connect Poisson generators if enabled
    use_poisson_generators=nest_pms["network"]["use_poisson_generators"]
    if use_poisson_generators:
        num_poisson_generators=nest_pms["network"]["num_poisson_generators"]
        poisson_rate=nest_pms["poisson"]["poisson_rate"]
        pgs = nest.Create("poisson_generator", num_poisson_generators, params={"rate": poisson_rate})
        poisson_weight=nest_pms["poisson"]["poisson_weight"]
        for i in range(num_exc_pop):
            if use_standard_exc_nest_neuron:
                nest.Connect(pgs, neurons[i],\
                             syn_spec={"weight": poisson_weight, "delay": 1.0})
        else:
            assert("MC neuron not yet implemented")           

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
            if use_standard_exc_nest_neuron:
                nest.Connect(dcgs[i], neurons[i],\
                    syn_spec={'weight': dc_exc_weight, 'delay': dc_exc_delay_ms})
            else:
                assert("MC neuron not yet implemented")    

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
        
    # Connect spike detectors
    for i in range(num_exc_pop):
        nest.Connect(neurons[i], spike_recorders[i])
    nest.Connect(inh_neurons, inh_spike_recorder)

    # Run the simulation
    nest.Simulate(sim_pms["stop_ms"])

    sim_completed=True
    return sim_completed, spike_recorders, inh_spike_recorder
