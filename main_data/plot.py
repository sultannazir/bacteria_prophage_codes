import numpy as np
import pandas as pd
import seaborn as sns
import math
import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec
import matplotlib.lines as mlines

R = 7 # Value of pi or RoC goes from 1 to R
N = 20 # Number of iterations goes from 1 to N

fig, axes = plt.subplots(1,3, figsize=(15, 6), tight_layout=True, dpi=300)
gs = GridSpec(1,3, figure=fig)

VDs=[3,1,-1]
VDtitles = [r"$V_{diff}=0.001$",r"$V_{diff}=0.1$","Virions well-mixed"]


for vidx in range(3):
    result = []
    bardata = []
    annot_values = []
    for r in range(0,R):
        E = 0
        for n in range(0,N):
            # if r%2==0:
            #     rr = str(int((r+8)/2))
            # else:
            #     rr = str((r+8)/2)
            dat = pd.read_csv('local_invasion/mix0timeseries_pi'+str(r+8)+'VD'+str(VDs[vidx])+'iter'+str(n+1)+'.dat', delimiter='\t', names=['time', 'U', 'Uc', 'L', 'Lc', 'Lp', 'Lcp', 'V', 'VT'])
            dat['fracC'] = (dat['Uc']+dat['Lc']+dat['Lcp'])/(dat['Uc']+dat['Lc']+dat['Lp']+2*dat['Lcp'])
            dat['fracP'] = (dat['Lp']+dat['Lcp'])/(dat['Uc']+dat['Lc']+dat['Lp']+2*dat['Lcp'])
            dat['gene'] = dat['Uc']+dat['Lc']+dat['Lp']+2*dat['Lcp']
            dat['genecell'] = dat['Uc']+dat['Lc']+dat['Lp']+2*dat['Lcp']

            frac = dat.fracP.iloc[-1]
            genepop = dat.gene.iloc[-1]
            if genepop==0:
                ext = "True"
            else:
                ext = "False"
                E += 1

            result.append([r,frac,ext,n+1])

        # Store the sum for annotation
        annot_values.append((r, E))

    data = pd.DataFrame(result, columns = ['pi', 'fraction', 'Gene extinction', 'replicate'])

    sns.violinplot(x='pi', y='fraction', data=data[data['Gene extinction']=="False"], cut=0,
               density_norm='width', inner=None, linewidth=1, color='#DDDDFF',
               saturation=1, ax = axes[vidx])

    g = sns.stripplot(x=data['pi'], y='fraction', data=data[data['Gene extinction']=="False"],
              jitter=True, linewidth=1, ax = axes[vidx], c='#2222AA', alpha=0.5, label=None)

    # Add annotations above each group
    for x, count in annot_values:
        axes[vidx].text(x, 1.03, str(count), ha='center', va='bottom', fontsize=12, color='black', weight="bold")

    axes[vidx].set_xlabel(r"Privatization (log-scaled), $log_{10}\pi$", fontsize=14)
    axes[vidx].set_ylabel("")
    axes[vidx].set_xlim(-1,R)
    axes[vidx].set_ylim(-0.05,1.15)
    axes[vidx].set_title(VDtitles[vidx], fontsize = 14)
    axes[vidx].grid()
    axes[vidx].set_xticks([0,1,2,3, 4, 5, 6])
    axes[vidx].set_xticklabels([-2.0, -2.25, -2.5, -2.75, -3.0, -3.25, -3.5])
    axes[vidx].invert_xaxis()
axes[0].set_ylabel("Fraction of prophage-association", fontsize=14)

# plt.show()

plt.savefig("local_invasion_frac_of_prophage_association_pi.png")
