def check_sim_recording_pms(sim_pms,recording_pms):
    assert(sim_pms["stop_ms"]>=0.0)
    assert(recording_pms["start_ms"]>=0.0)
    assert(recording_pms["stop_ms"]>=recording_pms["start_ms"])
    assert(recording_pms["start_ms"]<=sim_pms["stop_ms"])
    assert(recording_pms["duration_ms"]>=0.0)
    
def check_analysis_pms(crop_pms,recording_pms,analysis_pms):
    assert(crop_pms["stop_ms"]<=recording_pms["stop_ms"])
    assert(crop_pms["start_ms"]>=recording_pms["start_ms"])
    assert(crop_pms["duration_ms"] <= recording_pms["duration_ms"])
    assert(analysis_pms["spikes_sampling_window_ms"]<=crop_pms["duration_ms"])
    assert(analysis_pms["spectrogram_window_ms"]<=crop_pms["duration_ms"]/2)