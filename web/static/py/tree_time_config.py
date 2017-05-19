
concentration_params = [{"name": "[A]", "value":0.25}, {"name": "[C]", "value":0.25}, {"name": "[G]", "value":0.25}, {"name": "[T]", "value":0.25}]
available_gtrs =  [
            {"key":"jc", "value": "Jukes, Cantor 1969"},
            {"key":"k80", "value": 'Kimura 1980', 'params': [{"name":"kappa", "value":0.1}]},
            {"key":'f81', "value": 'Felsenstein 1981','params':concentration_params},
            {"key":'hky', "value": 'Hasegawa, Kishino, Yano 1985', 'params':concentration_params + [{"name":"kappa", "value":0.1}]},
            {"key":'t92', "value": 'Tamura 1992', 'params':[{"name":"[G+C]", "value":"0.3"}, {"name":"kappa", "value":0.1}]},
            {"key":'tn93', "value": 'Tamura, Nei 1993', 'params':concentration_params + [{"name":"kappa1", "value":0.1}, {"name":"kappa2", "value":0.5}]},
            {"key":'jtt', "value": 'Jones, Taylor, Thronton'},
        ]

treetime_webconfig = {
    # build tree using FastTree?
    "build_tree": True,
    # need to show confidence intervals later
    'do_marginal': True,
    'gtr': 'infer',
    'polytomies': True,
    'root': 'best',
    'slope': False,
    'slope_value' : 1e-3,
    'use_coalescent_prior':False,
    "coalescent_prior_value":0.01,
    'use_relaxed_clock':False,
    'relaxed_clock':
        {'slack' : 0.01,
         'coupling':0.005
        },
    'available_gtrs': available_gtrs
}



treeanc_webconfig = {
    # build tree using FastTree?
    "build_tree": True,
    # need to show confidence intervals later
    'do_marginal': True,
    'gtr': 'infer',
    'polytomies': True,
    'root': 'best',
    'available_gtrs': available_gtrs

}

if __name__ == "__main__":
    pass
