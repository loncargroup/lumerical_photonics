# Michael Haas (061523)
# msg_client.py
# Make a simplistic version of
from cavity_functions import make_cavity, simulate_cavity, gds_cavity

cp = {'a': 0.296 * 1e-6,
      'cX': 0.701,
      'cY': 1.7,
      'cW': 2,
      'cH': .5,
      'dA': 0.176,
      'dY': 0.2,
      'dX': 0.176,
      'dA_tap': 0,
      'dY_tap': 1,
      'dX_tap': 0,
      'nleft': 6,
      'nright': 6,
      'ndef': 6,
      'ntaper': 10,
      'hole_type': 'rib',
      'shift': -1,
      'source frequency': 406.7e12,
      # 'ref_index': 2.4028,
      'effective_index': 1.6,
      'resonance_wavelength': 0.737,
      'cav_defect_type': 'cubic',
      'tap_defect_type': 'cubic',
      'min_dim': 5e-8
      }

cavity = make_cavity(cp)  # This works for the cavity, store the hole list as well
# cav_results = simulate_cavity(cavity, cp)
gds = gds_cavity(cavity, cp)
