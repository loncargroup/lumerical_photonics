from wvgsolver.geometry import MeshRegion
from wvgsolver import Vec3
from wvgsolver.utils import BBox
import numpy as np


def simulate_cavity(cavity, cp):
    source_frequency = cp['source frequency']
    man_mesh = MeshRegion(BBox(Vec3(0), Vec3(12e-6, 0.7e-6, 0.4e-6)), 12e-9, dy=None, dz=None)

    # There are a few different simulation conditions, let's keep this one but see if we need to change later
    r1 = cavity.simulate("resonance", target_freq=source_frequency, source_pulselength=9e-15, analyze_fspan=10e12,
                         analyze_time=500e-15, mesh_regions=[man_mesh],
                         sim_size=Vec3(1.25, 3, 10))  # before 600 analyze time and 60 source pulse length
    # r1 = cavity.simulate("resonance", target_freq=target_frequency, mesh_regions=[man_mesh],
    #                     sim_size=Vec3(2, 3, 8))

    # r1 = cavity.simulate("resonance", target_freq=source_frequency)
    if cp['plot_E']:
        r1["xyprofile"].show()
        r1["yzprofile"].show()

    c = 2.99792458e17  # This is the speed of light in nm/s. As such, in order to get wavelength, it's just dividing c/freq
    qleft = r1["qxmin"]
    qright = r1["qxmax"]
    qy = 1 / (2 / r1["qymax"])
    qz = 1 / (1 / r1["qzmin"] + 1 / r1["qzmax"])
    freq = r1["freq"]
    vmode = r1["vmode"]

    qperp = 1 / (1 / qy + 1 / qz)
    qscat = 1 / (1 / qperp + 1 / qleft)
    qwvg = qright

    wavelen_pen = np.exp(-((source_frequency - freq) / 4e12) ** 2)
    rerun_thresh = 0.90

    if wavelen_pen < rerun_thresh:
        # shift source frequency to cavity resonance and rerun simulation.
        # (this should help avoid non cavities with artificially low mode volumes)
        print('Source frequency =' + str(source_frequency))
        print('Simmed frequency =' + str(freq))
        print('Wavelen_penalty =' + str(wavelen_pen))
        print('---n')
        cp['source frequency'] = freq
        witness_rerun = simulate_cavity(cavity, cp)
        return witness_rerun

    cavity_results = {
        # Here are the values that we'll optimize over
        'qscat': qscat,
        'qwvg': qwvg,
        'qleft': qleft,
        'vmode': vmode,
        # This is to keep track of wavelengths
        'freq': freq,
        'wl': c / freq,
    }

    print('---')
    return cavity_results
