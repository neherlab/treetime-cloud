config_dic = {
    "use_mu":False,
    "mu": 1e-3,
    "gtr": "infer", # possible values: infer, Jukes-Cantor
    "resolve_poly": False,
    "build_tree": False,
    "model_coalescent":False,
    "Tc": 0.01, # Coalescent model Tc=0 switches the option off, otherwise use t
                # his value as coalescence timescale
    "reroot": "best",
    "do_relaxed_clock" : False,
    "relax_clock": {"slack": 0, "coupling": 0}, # slack=0 && coupling==0 switch the
                                                # option to compute relaxed clock OFF.
                                                # otherwise, their values are used to
                                                # relax the molecular clock
    "reuse_branch_len": True
    }

if __name__ == "__main__":
    pass