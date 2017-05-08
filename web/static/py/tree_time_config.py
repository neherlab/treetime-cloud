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
    'available_gtrs':
        [
            {"key":"jc", "value": "Jukes, Cantor 1969"},
            {"key":"k80", "value": 'Kimura 1980', 'params': ['kappa']},
            {"key":'f81', "value": 'Felsenstein 1981','params':['Ca', 'Cc', 'Cg', 'Ct']},
            {"key":'hky', "value": 'Hasegawa, Kishino, Yano 1985', 'params':['Ca', 'Cc', 'Cg', 'Ct', 'kappa']},
            {"key":'t92', "value": 'Tamura 1992', 'params':['Cgc', 'kappa']},
            {"key":'tn93', "value": 'Tamura, Nei 1993', 'params':['Ca', 'Cc', 'Cg', 'Ct', 'kappa1', 'kappa2']},
            {"key":'jtt', "value": 'Jones, Taylor, Thronton'},
        ]
}

if __name__ == "__main__":
    pass
