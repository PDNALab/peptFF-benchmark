# GROMACS Simulation Workflow

This project provides a structured pipeline to run molecular dynamics simulations in **GROMACS 2019.2**, using various force fields and corresponding water models. 
The workflow includes system preparation and processing trajectory.

---

## Software Dependencies

The following modules or software versions were loaded to run GROMACS simulations:

- `cuda/10.0.130`
- `intel/2018.1.163`
- `openmpi/4.0.0`
- `gromacs/2019.2`

---

## Supported Force Fields & Water Models

| Force Field        | Directory/Water Model                  |
|--------------------|----------------------------------------|
| a99SBdisp          | `a99SBdisp.ff/a99SBdisp_water.gro`     |
| charmm36IDPSFF     | `tip3p.gro`                            |
| des-amber          | `tip4p.gro`                            |
| des-amber-SF1.0    | `tip4p.gro`                            |
| OPLSIDPSFF         | `tip4p.gro`                            |

Ensure the relevant force field directory is present in your working directory.

---

## Input Files

- `peptide.pdb`: Initial structure file with hydrogens (from **xleap**).
- `peptide_reduced.pdb`: Same structure without hydrogens.

---

## MDP Files Required

| Filename     | Description                                  |
|--------------|----------------------------------------------|
| `ions.mdp`   | Add ions to system                           |
| `min1.mdp`   | Energy minimization (Steepest Descent)       |
| `min2.mdp`   | Energy minimization (Conjugate Gradient)     |
| `nvt.mdp`    | NVT equilibration (10 ns)                    |
| `npt.mdp`    | NPT equilibration (20 ns)                    |
| `md.mdp`     | Production run (200 ns, NPT conditions)      |

Use the provided script (`simulation.sh`) to execute the pipeline.

---

## Simulation Workflow (Example: a99SBdisp.ff)

```bash
# Step 1: Generate topology
gmx pdb2gmx -f peptide_reduced.pdb -o peptide_processed.gro

# Step 2: Define simulation box
gmx editconf -f peptide_processed.gro -o peptide_newbox.gro -c -d 1.0 -bt octahedron

# Step 3: Solvate system
gmx solvate -cp peptide_newbox.gro -cs a99SBdisp.ff/a99SBdisp_water.gro -o peptide_solv.gro -p topol.top

# Step 4: Prepare for ion addition
gmx grompp -f ions.mdp -c peptide_solv.gro -p topol.top -o ions.tpr

# Step 5: Add ions
gmx genion -s ions.tpr -o peptide_solv_ions.gro -p topol.top -pname NA -nname CL -neutral -conc 0.1

Note: check topol.top at every step to avoid errors.



