from cavity_functions import make_cavity, simulate_cavity
from scipy import optimize
import pandas as pd


def fitness(opt_params):
    global dat_loc
    global df
    cp = {
        'a': opt_params[0] * 1e-6,
        'w0': opt_params[1] * 1e-6,
        'h0': opt_params[2] * 1e-6,
        'cX': opt_params[3],
        'cY': opt_params[4],
        'dA': opt_params[5],
    }

    constants = {
        'dX': cp['dA'],
        'dY': 0,
        'nleft': 10,
        'ndef': 6,
        'nright': 3,
        'ntaper': 3,
        'dA_tap': 0,
        'dY_tap': .95,
        'dX_tap': 0,
        'hole_type': 'rib',
        'shift': -1,
        'source frequency': 406.7e12,
        'target frequency': 406.7e12,
        'cav defect type': 'cubic',
        'plot_E': False,
        'min dim': .05e-6
    }

    cp.update(constants)
    cavity = make_cavity(cp)
    if cavity is None:
        return 0
    cav_results = simulate_cavity(cavity, cp)

    qscat = cav_results['qscat']
    vmode = cav_results['vmode']
    fitness_val = -qscat / vmode
    data_dict = {'a': cp['a'],
                 'w0': cp['w0'],
                 'h0': cp['h0'],
                 'cX': cp['cX'],
                 'cY': cp['cY'],
                 'dA': cp['dA'],
                 'Waveguide Q': cav_results['qwvg'],
                 'Scattering Q': cav_results['qscat'],
                 'Mode Volume': cav_results['vmode'],
                 'Wavelength': cav_results['wl'],
                 'Fitness Value': -fitness_val
                 }
    df = df.loc[len(df)] = data_dict
    df.to_csv(dat_loc)
    return fitness_val


# [a, w0, h0, cX, cY, dA]
starting_values = [0.292, .6, .15, .7, 1.7, .17]
lb = [.25, .3, .1, .5, 1, 0]
ub = [.4, .8, .3, .9, 2, .3]

bounds = optimize.Bounds(lb, ub)
directory = 'G:/My Drive/Grad School/Simulations/lum_optimization'
filename = '062823_riboptimization_a-w0-h0-cX-cY-dA.csv'

dat_loc = directory + '/' + filename
column_data = ['a', 'w0', 'h0', 'cX', 'cY', 'dA', 'Waveguide Q',
               'Scattering Q', 'Mode Volume', 'Wavelength', 'Fitness value']
df = pd.DataFrame(columns=column_data)
df.to_csv(dat_loc)
# res = optimize.minimize(fun=fitness, args = df, x0=starting_values, bounds=bounds)
res = fitness(starting_values)
