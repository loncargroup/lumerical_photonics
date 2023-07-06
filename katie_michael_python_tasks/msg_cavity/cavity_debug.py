from cavity_functions import make_cavity, simulate_cavity

cp = {  # Actual dimensions
    'a': 0.2965 * 1e-6,
    'w0': .6 * 1e-6,
    'h0': .15 * 1e-6,
    # Scalars of lattice constant
    'cX': 0.701,
    'cY': 1.7,
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
    'source frequency': 406.7e12,
    'cav_defect_type': 'cubic',
    'plot_E': True
}

cavity = make_cavity(cp)
cav_results = simulate_cavity(cavity, cp)
print(cav_results)
