from make_functions import make_cavity
from simulate_functions import simulate_cavity
from scipy.optimize import minimize, LinearConstraint
import numpy as np


def fitness(params):
    global source_frequency
    # print(params)  # <-- you'll see that params is a NumPy array
    aIn, w0In, h0In, hxIn, hyIn = params  # <-- for readability you may wish to assign names to the component variables
    print(params)
    cp = {  # Actual dimensions
        'a': aIn * 1e-6,
        'w0': w0In * 1e-6,
        'h0': h0In * 1e-6,
        # Scalars of lattice constant
        'hx': hxIn * 1e-6,
        'hy': hyIn * 1e-6,
        'dA': 0.176,
        'dX': 0.176,
        'dY': 0,
        'dA_tap': 0,
        'dY_tap': .95,
        'dX_tap': 0,
        # Mirror cell definitions
        'nleft': 10,
        'ndef': 6,
        'nright': 3,
        'ntaper': 3,
        # Cavity specifics
        'hole_type': 'rib',
        'shift': -1,
        'source frequency': source_frequency,
        'cav defect type': 'cubic',
        'min dim': .05 * 1e-6,
        'plot_E': False
    }

    cavity, flag = make_cavity(cp)
    if flag:
        cav_results = simulate_cavity(cavity, cp)
        qsc = cav_results['qscat']
        vmode = cav_results['qwvg']
        print(cav_results)
        source_frequency = cav_results['freq']
        return -qsc / vmode
    else:
        return 0


source_frequency = 406.7 * 1e12
initial_guess = [.295, .6, .15, .2, .5]
min_dim = .05
# linear_constraint = LinearConstraint([[1, 0, 0, -1, 0]
#                                      , [0, 0, 0, 1, 0]
#                                      , [0, 0, 0, 0, 1]
#                                      , [0, 1, 0, 0, -1]
#                                      , [0, 0, 1, 0, 0]], np.full((1, 5), min_dim)[0])

result = minimize(fitness, initial_guess, method='Nelder-Mead')
