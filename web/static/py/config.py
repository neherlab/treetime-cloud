config_dic = {
    "mu": 0.0,
    "gtr": "infer", # possible values: infer, Jukes-Cantor
    "resolve_poly": True,
    "build_tree": False,
    "Tc": 0.01, # Coalescent model Tc=0 switches the option off, otherwise use t
                # his value as coalescence timescale
    "reroot": "best",
    "relax_clock": {"slack": 0, "coupling": 0}, # slack=0 && coupling==0 switch the
                                                # option to compute relaxed clock OFF.
                                                # otherwise, their values are used to
                                                # relax the molecular clock
    "reuse_branch_len": True
    }


if __name__ == "__main__":
    pass