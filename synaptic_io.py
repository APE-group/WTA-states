import pickle

def load_intra_assembly_syn(file_name, verbose):
    if(verbose):
        print("in load_intra_assembly_syns:")
        print("reading from file:", file_name)
    with open(file_name, "rb") as f:
            array_of_dicts = pickle.load(f)
    return array_of_dicts
