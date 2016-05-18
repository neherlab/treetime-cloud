config_dic = {
    #  set the configuration parameters
    "do_build_tree":True,
    "reuse_branch_len":True,  # Do we trust the branch lengths of the tree given? 
    "infer_gtr":True,  # Infer GTR model from the data given (Use default if False)?
    "gtr":"jukes_cantor",
    "reroot": True,  # Find the root node which maximizes the molecular clock correlation?
    "use_mu":False,  # Should we use the mutation rate (if False, it will be inferred from the molecular clock)? 
    "mu":None,  #  #mutations per years per position
    "resolve_poly":True,  # Try to resolve multiple mergers
    "coalescent": False,  # Run the coalescent model (takes care about polytomies resolution)
    "Tc":0.00,  #  coalescent distance (%difference )
    "relax_mu":False,  #  Relax mutation rate? 
    "slack":0.01,  # How much we allow for variation between parents and children 
    "coupling":0.002
}

if __name__ == "__main__":
    pass