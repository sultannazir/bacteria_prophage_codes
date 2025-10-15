import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
from matplotlib.gridspec import GridSpec

N = 20 # Number of iterations goes from 1 to N

fig, axes = plt.subplots(7,3, figsize=(10, 13), tight_layout=True, dpi=300)
gs = GridSpec(7,3, figure=fig)

VDs=[3, 1, -1]
VDtitles = [r"$V_{diff}=0.001$",r"$V_{diff}=0.1$","Virions well-mixed"]



for pi in range(7):
    for vidx in range(3):
        for n in range(0,N):
            if pi%2==0:
                rr = str(int((pi+8)/2))
            else:
                rr = str((pi+8)/2)
            dat = pd.read_csv('data/mix0pi'+rr+'VD'+str(VDs[vidx])+'iter'+str(n+1)+'IwGR0timeseries.dat', delimiter='\t', names=['time', 'U', 'Uc', 'L', 'Lc', 'Lp', 'Lcp', 'V', 'VT'])
            dat['fracP'] = (dat['Lp']+dat['Lcp'])/(dat['Uc']+dat['Lc']+dat['Lp']+2*dat['Lcp'])
            axes[pi][vidx].plot(dat['fracP'])
            axes[pi][vidx].set_ylim(-0.1,1.1)
    aa = axes[pi][-1].twinx()
    aa.set_yticks([])
    aa.set_ylabel(r"$\log_{10} \pi = $"+str(-(pi+8)/4), rotation=0, labelpad = 40)

for vidx in range(3):
    axes[0][vidx].set_title(VDtitles[vidx])
    axes[-1][vidx].set_xlabel("Time "+r"$/10^{3}$")
axes[3][0].set_ylabel("Fraction of prophage-association")

#plt.show()
plt.savefig("timeseries.png")
