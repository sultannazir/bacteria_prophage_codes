import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec
import pandas as pd
import numpy as np

pi = 4
VD = 3

grid_time = 125000

genotype_to_colour = {-1: np.array([0, 0, 0]),
                      'U': np.array([255, 251, 158]),
                      'Uc': np.array([111, 46, 156]),
                      'L': np.array([247, 172, 70]),
                      'Lc': np.array([66, 0, 100]),
                      'Lp': np.array([206, 80, 50]),
                      'Lcp': np.array([20, 0, 68])}

timeseries = pd.read_csv('timeseries_data/pi'+str(pi)+'VD'+str(VD)+'/timeseries.dat', delimiter='\t', names=['time', 'U', 'Uc', 'L', 'Lc', 'Lp', 'Lcp', 'V', 'VT'])

grid = pd.read_csv('grids/pi'+str(pi)+'VD'+str(VD)+'/grid_'+str(grid_time)+'.dat', delimiter='\t', names = ['x', 'y', 'genotype'])
grid = grid.pivot(index='x', columns='y', values='genotype')

# Replace '-1' string with actual -1 key
grid = grid.replace('-1', -1)

# Map each element to its RGB array using .apply + .map
mapped = grid.apply(lambda col: col.map(lambda x: genotype_to_colour[x]))

fig = plt.figure(constrained_layout=True, figsize=(8,8), dpi=300)
gs = GridSpec(3,2, figure=fig)

ax1 = fig.add_subplot(gs[0:2,:])
ax2 = fig.add_subplot(gs[2,:])

ax1.imshow(
    np.array(
        [[np.array(i, dtype=np.uint8)for i in j] for j in mapped.values],
        dtype=np.uint8,
    )
)
ax1.axis('off')

ax2.plot(timeseries['time'], timeseries['U'], color='#FFFB9E', linewidth=3)
ax2.plot(timeseries['time'], timeseries['Uc'], color='#6F2E9C', linewidth=3)
ax2.plot(timeseries['time'], timeseries['L'], color='#F7AC46', linewidth=3)
ax2.plot(timeseries['time'], timeseries['Lc'], color='#420064', linewidth=3)
ax2.plot(timeseries['time'], timeseries['Lp'], color='#CE5035', linewidth=3)
ax2.plot(timeseries['time'], timeseries['Lcp'], color='#140044', linewidth=3)
ax2.grid()
ax2.spines[['right', 'top']].set_visible(False)
ax2.set_xlim(0,250000)
ax2.set_ylim(0,15000)
ax2.set_xlabel('Time', fontsize=17)
ax2.set_ylabel('Abundance', fontsize=17)
ax2.set_yticks([0, 5000, 10000])
ax2.axvline(grid_time, c='black', linewidth=5, alpha=0.5)
#plt.show()
plt.savefig('grids/figures/grid_pi4VD3.png')
