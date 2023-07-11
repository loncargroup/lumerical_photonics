import pandas
import matplotlib.pyplot as plt

directory = 'G:/Shared drives/SiV-OMC/gds'
fname = '070623_ribparametersweep_a-cY_simmed'

file = directory + '/' + fname + '.csv'

df = pandas.read_csv(file)
c = 299792.458 * 1e3  # nm * THz
wl = [c / i for i in df['Frequency']]

fig = plt.figure()
ax = fig.add_subplot(projection='3d')

# for i in range(95):
#     ax.scatter(df['a'][i], df['cY'][i], df['Scattering Q'][i])
# plt.show()

for i in range(95):
    ax.scatter(df['a'][i], df['cY'][i], wl[i])
plt.show()
