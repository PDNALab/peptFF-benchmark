import MDAnalysis as mda
from MDAnalysis.analysis.dihedrals import Ramachandran
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec

plt.switch_backend('agg')

def extract_angles(topology, trajectory):
    u = mda.Universe(topology, trajectory)
    r = u.select_atoms('backbone')
    R = Ramachandran(r).run()
    phi_angles = R.angles[:, :, 0].flatten()
    psi_angles = R.angles[:, :, 1].flatten()
    return phi_angles, psi_angles

# Extract angles
angle_data = [
    extract_angles("ff19SB/stripped.peptide_sol.prmtop", "ff19SB/stripped_md.nc"),
    extract_angles("ff99SB/stripped.peptide_sol.prmtop", "ff99SB/stripped_md.nc"),
    extract_angles("ff99IDPs/stripped.peptide_sol.prmtop", "ff99IDPs/stripped_md.nc"),
    extract_angles("ff14IDPs/stripped.peptide_sol.prmtop", "ff14IDPs/stripped_md.nc"),
    extract_angles("ff14IDPSFF/stripped.peptide_sol.prmtop", "ff14IDPSFF/stripped_md.nc"),
    extract_angles("OPLSIDPSFF/production/md_stripped.pdb", "OPLSIDPSFF/production/md_stripped.xtc"),
    extract_angles("a99SBdisp/production/md_stripped.pdb", "a99SBdisp/production/md_stripped.xtc"),
    extract_angles("des-amber/production/md_stripped.pdb", "des-amber/production/md_stripped.xtc"),
    extract_angles("des-amber-SF1.0/production/md_stripped.pdb", "des-amber-SF1.0/production/md_stripped.xtc"),
    extract_angles("charmm36m/stripped.peptide_sol.pdb", "charmm36m/stripped_md.nc"),
    extract_angles("charmm36IDPSFF/production/md_stripped.pdb", "charmm36IDPSFF/production/md_stripped.xtc"),
]

# Create free energy landscapes
def create_free_energy(phi, psi):
    hist, xedges, yedges = np.histogram2d(phi, psi, bins=360, range=[[-180, 180], [-180, 180]], density=True)
    with np.errstate(divide='ignore'):
        F = -np.log(hist)
        F -= np.min(F)
    return np.ma.masked_where(np.isinf(F), F)

free_energies = [create_free_energy(phi, psi) for phi, psi in angle_data]

# Forcefield labels (only 11)
forcefields = [
    'ff19SB', 'ff99SB', 'ff99IDPs', 'ff14IDPs', 'ff14IDPSFF',
    'OPLSIDPSFF', 'a99SBdisp', 'des-amber', 'des-amber-SF1.0',
    'charmm36m', 'charmm36IDPSFF'
]

# Set up figure with 12 subplots (4x3)
fig = plt.figure(figsize=(15, 15))
gs = GridSpec(4, 3, hspace=0, wspace=0.2)
axs = [fig.add_subplot(gs[i]) for i in range(12)]

# Plot only the first 11 with data
for i in range(11):
    im = axs[i].imshow(free_energies[i].T, origin='lower', extent=[-180, 180, -180, 180], cmap='RdBu', aspect='auto')
    axs[i].set_ylabel(forcefields[i], fontsize=20)

# Turn off the 12th subplot
axs[11].axis('off')

# Set x and y ticks
for ax in axs[:11]:
    ax.set_xticks([-180, 0, 180])
    ax.set_yticks([-180, 0, 180])
    ax.set_xticklabels([-180, 0, 180], fontsize=16)
    ax.set_yticklabels([-180, 0, 180], fontsize=16)

# Remove x-tick labels from non-bottom rows
#for ax in axs[:-3]:
#    ax.set_xticklabels([])

# Remove x-axis tick labels from all except last row + subplot 5
for i, ax in enumerate(axs[:11]):
    if i not in [8, 9, 10]:  # Keep labels only for these
        ax.set_xticklabels([])

# Set x-axis labels for last row and subplot 5
#for i in [8, 9, 10]:
#    axs[i].set_xlabel("Phi (°)", fontsize=20)

for i, ax in enumerate(axs[:11]):
    if i % 3 != 0:  # Not in first column
        ax.set_yticklabels([])

for i in range(11):
    axs[i].set_ylabel(forcefields[i], fontsize=20)

# Colorbar
cbar_ax = fig.add_axes([0.92, 0.15, 0.02, 0.7])
cbar = fig.colorbar(axs[4].images[0], cax=cbar_ax, label='Free Energy (kT)')
cbar.set_label('Free Energy (kT)', fontsize=20)
cbar.ax.tick_params(labelsize=20)

# Save
plt.tight_layout(rect=[0.12, 0.06, 0.9, 1])
#plt.tight_layout(rect=[0.06, 0.06, 0.9, 1])
plt.savefig('mda_Ramachandran_all.png', transparent=True, dpi=600)
plt.show()

