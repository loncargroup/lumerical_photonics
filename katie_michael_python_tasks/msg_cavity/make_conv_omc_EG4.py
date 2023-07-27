# # The goal of this is to put our conventional omc design onto the new egrade pattern. We have easy code to generate
# the egrade, now it's just a matter of putting our old omc design onto the layout.
from dev_layout import dev_map
from old_omc_functions import deviceInverse2
import numpy as np

layout = {'layout_index': 0, 'aper_index': 1, 'layer': 4, 'directory': 'G:/Shared drives/SiV-OMC/gds',
          'filename': '65nm_EG4_omc'}

column_data = ["Device x position", "Device y position", "a", "cX"]

dat_loc = layout['directory'] + '/' + layout['filename']
gds_name = dat_loc + '.gds'
df_name = dat_loc + '.csv'
E, pos_list = dev_map(layout)
# Let's convert over the old omc design.
i = 0
num_d = 9
num_w = 6
dArr = np.linspace(.085,.105,num_d)
dMArr = [i + .02 for i in dArr]
WmList = np.linspace(.235,.267, num_w)
NrightList = [7,10]

for (d0, dM) in zip(dArr,dMArr):
    for Wm in WmList:
        for Nright in NrightList:
            (x, y) = pos_list[i]
            a = deviceInverse2(x, y, d0, dM, Nright, Wm, i, np.pi / 4)
            i+=1
            E.add(a)

E.write_gds(gds_name)
