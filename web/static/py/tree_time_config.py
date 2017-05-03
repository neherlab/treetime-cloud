treetime_webconfig = {
    # build tree using FastTree?
    "build_tree": False,
    # need to show confidence intervals later
    'do_marginal': True,
    'gtr': 'infer',
    'polytomies': True,
    'root': 'best',
    'slope': None,
    'relaxed_clock':False,
    'coalescent':None,
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
