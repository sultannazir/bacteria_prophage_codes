import numpy as np
import pandas as pd
import seaborn as sns
import math
import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec
import matplotlib.lines as mlines

R = 6 # Value of pi or RoC goes from 1 to R
N = 20 # Number of iterations goes from 1 to N

fig, axes = plt.subplots(1,3, figsize=(13, 6), tight_layout=True, dpi=300)
gs = GridSpec(1,3, figure=fig)

VDs=[-1,1,3]
VDtitles = ["Virions well-mixed",r"$V_{diff}=0.1$",r"$V_{diff}=0.001$"]


for vidx in range(3):
    result = []
    bardata = []
    for r in range(0,R):
        nFnE = 0
        FnE = 0
        nFE = 0
        FE = 0
        for n in range(0,N):
            dat = pd.read_csv('pi_varies_selection/cells_stationary/R1pi'+str(r+1)+'iter'+str(n+1)+'base_VD'+str(VDs[vidx])+'.dat', delimiter='\t', names=['time', 'U', 'Uc', 'L', 'Lc', 'Lp', 'Lcp', 'V', 'VT'])
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

            result.append([r,frac,ext,n+1])

    data = pd.DataFrame(result, columns = ['pi', 'fraction', 'Gene extinction', 'replicate'])

    sns.violinplot(x='pi', y='fraction', data=data[data['Gene extinction']=="False"], cut=0,
               density_norm='width', inner=None, linewidth=1, color='#DDDDFF',
               saturation=1, ax = axes[vidx])

    g = sns.stripplot(x=data['pi'], y='fraction', data=data[data['Gene extinction']=="False"],
              jitter=True, linewidth=1, ax = axes[vidx], c='#2222AA', alpha=0.5, label=None)

    axes[vidx].set_xlabel(r"Privatization (log-scaled), $log_{10}\pi$", fontsize=14)
    axes[vidx].set_ylabel("")
    axes[vidx].set_xlim(-1,R+1)
    axes[vidx].set_ylim(-0.05,1.05)
    axes[vidx].set_title(VDtitles[vidx], fontsize = 14)
    axes[vidx].grid()
    axes[vidx].set_xticks([0,1,2,3,4,5,6])
    axes[vidx].set_xticklabels([-0.5, -1.0, -1.5, -2.0, -2.5, -3.0, -3.5])
    axes[vidx].invert_xaxis()
axes[0].set_ylabel("Fraction of prophage-association", fontsize=14)

plt.savefig("frac_of_prophage_association_pi.png")
