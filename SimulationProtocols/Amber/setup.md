# AMBER Simulation Workflow

This project provides a structured pipeline to run molecular dynamics simulations in **AMBER 2020**, using various force fields and corresponding water models. 
The workflow includes system preparation and processing trajectory.

## Software Dependencies

The following modules or software versions were loaded to run GROMACS simulations:

- `python3`
- `cuda/11.1.0`
- `nvhpc/20.11`
- `openmpi/4.0.5`
- `amber/20`


## Supported Force Fields & Water Models & CMAPs

| Force Field        | Water Mod                   | CMAP              |
|--------------------|-----------------------------|-------------------|
| ff19SB             | `OPC`                       |    `-`            |
| ff99SB             | `OPC`                       |    `-`            |
| ff99IDPs           | `OPC`                       | `ff99IDPs.para`   |
| ff14IDPs           | `TIP3P`                     | `ff14IDPs.para`   |
| ff14IDPsFF         | `TIP3P`                     | `ff14IDPsFF.para` |
| Charmm36m          | `TIP3P`                     |    `-`            |

To add the CMAP parameters corrections, follow the steps from: **[chaohao2010/ADD-CMAP](https://github.com/chaohao2010/ADD-CMAP)**.

For Charmm36m, topologies and coordinates were generated using: **[CHARMM-GUI](https://www.charmm-gui.org/)**

## Input Files

- `peptide.pdb`: Initial structure file (from **xleap**).
- `ADD-CMAP.py`: Add CMAP corrections


## MDP Files Required

| Filename     | Description                                           |
|--------------|-------------------------------------------------------|
| `tleap.in`   | solvate and dad ions to system                        |
| `min1-5.in`  | Energy minimization (with dereasing water restraints) |
| `min6.in`    | Energy minimization (no water restraints)             |
| `mdt.in`     | NVT equilibration (10 ns)                             |
| `npt1.in`    | NPT equilibration (10 ns)                             |
| `npt2.in`    | NPT equilibration (10 ns)                             |
| `md.in`      | Production run (200 ns, NPT conditions)               |

Use the provided script (`amber_sim.sh`) to execute the pipeline.

## Simulation Workflow (Example: ff14IDPs)

```bash
# Step 1: Generate topology & atomic coordinates
tleap -f tleap.in

# Step 2: Randomize ions positions
cpptraj -p peptide_sol.prmtop -i cpptraj_randomizeions.in

# Step 3: Add CMAP corrections
# Follow steps from https://github.com/chaohao2010/ADD-CMAP
python3 ADD-CMAP.py -p peptide_sol.prmtop -c ff14IDPs.para -o amber_CMAP.prmtop -s

# Step 4: Generate Hydrogen mass repartitioned topology
parmed -p peptide_CMAP.prmtop -i parm.in


Use the provided script (`amber_sim.sh`) to execute the run minimization, heating, equilibration and production simulations.

```
## Post-Processing Trajectories

```bash
# Merge all trajectories & remove water molecules and ions
cpptraj -i autoimage.in

```

All the analysis were performed using cpptraj and MDAnalysis tool. 


