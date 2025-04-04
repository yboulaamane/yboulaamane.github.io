---
title: 'An In-Depth Guide to MD Simulation and Analysis of a Protein–Ligand Complex'
date: 2025-04-04
permalink: /blog/2025/04/an-in-depth-guide-to-md-simulation-and-analysis-of-a-protein–ligand-complex
excerpt_separator: <!--more-->
toc: true
tags:
  - computational chemistry
  - drug discovery
---

Molecular Dynamics (MD) simulations allow computational chemists to explore the dynamic behaviour of protein–ligand complexes beyond static structures.
<!--more-->

# Introduction: Why Use MD for Protein–Ligand Systems?

A crystal or docked structure provides a single snapshot, but MD **bridges the gap between static protein–ligand structures and the dynamical changes** that govern biological function​ (1). By integrating Newtonian mechanics, MD generates an ensemble of conformations over time, revealing how a ligand binds, how a protein’s binding pocket adapts, and what stable or transient interactions form. In drug design and mechanistic studies, MD can help evaluate **binding stability, conformational flexibility, and key interactions** (hydrogen bonds, hydrophobic contacts, etc.) throughout a simulation. In short, MD simulations serve as a “computational microscope,” capturing protein and ligand motions at atomic resolution over nanoseconds to microseconds, thereby providing insights into **stability, conformational changes, and energetics** that are difficult to obtain experimentally.

**Why 100 ns?** A 100-nanosecond production run is a moderate timescale that is often sufficient for a protein–ligand complex to relax and sample local conformational changes (like side-chain reorientations or loop motions). While some slow events (large domain movements or complete unbinding) might require longer simulations, 100 ns strikes a practical balance between computational cost and capturing the complex’s equilibrium behaviour. In this example pipeline, we outline how to set up and analyse a 100 ns explicit-solvent MD simulation of a protein–ligand complex. The target audience is intermediate to advanced computational chemists, so we assume basic familiarity with MD concepts and GROMACS, focusing on **technical specifics, best practices, and scientific interpretation** of each analysis step.

# Simulation Setup: CHARMM36 Force Field, TIP3P Water, and System Preparation

**Force field and water model:** We use the **CHARMM36 all-atom force field** for the protein (and a compatible CHARMM-derived ligand parameter set), a well-validated choice for biomolecular simulations​ (2). CHARMM36 is known for accurately modelling protein secondary structure and side-chain interactions. It is designed to be used with the **TIP3P water model**, which is a three-site rigid water model. In fact, the CHARMM force fields employ a slightly modified TIP3P model (with Lennard-Jones parameters on hydrogens) to better reproduce hydration effects​ (3). Using TIP3P water ensures compatibility with CHARMM36 (since the force field was parametrized and tested with this water model) and provides a reasonable balance of accuracy and computational efficiency.

**Initial coordinates:** Typically, one starts from an experimental protein structure (e.g. a PDB file) with a bound ligand or a docked complex. Prepare separate coordinate files for the protein and ligand if needed. **Parameterization of the ligand** is crucial – here, since we use CHARMM36, one would obtain ligand parameters via the CHARMM General Force Field (CGenFF) (4) or a CHARMM-compatible source (5). Ensure that the ligand’s topology (charges, atom types, bonds, etc.) is generated and merged with the protein topology in GROMACS format​ (5).

**System preparation steps:** Setting up a protein–ligand simulation involves several standard steps​ (6):

- _Add missing atoms:_ Use tools (like gmx pdb2gmx in GROMACS or external programs like PDBFixer) to add any missing hydrogens or heavy atoms in the protein structure (e.g., missing loops or side chains), and assign appropriate protonation states for ionizable residues (consider pH). This yields a complete protein structure ready for the force field.

```
pdbfixer protein.pdb --output=protein_fixed.pdb
```

- _Apply force fields:_ Here, we save a clean copy of the protein structure without the ligand to generate the GROMACS topology using gmx pdb2gmx with the CHARMM36 force field. The ligand topology, generated separately via CGenFF, is then incorporated. At this stage, the full system topology—containing all necessary force field parameters for both protein and ligand—is assembled and ready for simulation.

```
grep UNK protein_fixed.pdb > unk.pdb
grep -v "UNK" protein_fixed.pdb > clean.pdb
gmx pdb2gmx -f clean.pdb -o processed.gro -ter -ignh
```

- _Generate ligand topology:_ To generate a topology for a small-molecule ligand compatible with the CHARMM36 force field, begin by saving the ligand in **.mol2** format with correct atom types and explicit hydrogens using tools like **Avogadro** or **OpenBabel**. Next, submit the ligand structure to the **CGenFF (ParamChem) server**. Upon successful submission, the server returns a **.str file** (CHARMM stream file) containing the ligand’s topology and force field parameters.

```
obabel unk.pdb -O unk.mol2 --addh --gen3d
sed -i '2s/.\*/UNK/; s/UNK900/UNK/g' unk.mol2
perl sort_mol2_bonds.pl unk.mol2 unk_fix.mol2
~ /silcsbio.2024.1/cgenff/cgenff unk_fix.mol2 -f unk.str (or via the web server)
python cgenff_charmm2gmx.py UNK unk_fix.mol2 unk.str charmm36-jul2022.ff
gmx editconf -f unk_ini.pdb -o unk.gro
```

- Merge the protein and ligand: To prepare the simulation system, merge the processed protein (processed.gro) and ligand (unk.gro) into a single GROMACS structure file:

```
cp processed.gro complex.gro && sed -i '$d' complex.gro && tail -n +3 unk.gro >> complex.gro && echo " " >> complex.gro
num1=$(sed -n '2p' complex.gro)
num2=$(sed -n '2p' unk.gro)
sum=$((num1 + num2))
sed -i "2s/.\*/ $sum/" complex.gro
```

- Update topology file: Update topol.top to include the ligand topology and parameters:

```
sed -i '/; Include Position restraint file/{:a;N;/#endif/!ba;s/\\n#endif\\n/\\n#endif/;s/#endif/#endif\\n\\n; Include ligand topology\\n#include "unk.itp"/}' topol.top
sed -i '/; Include forcefield parameters/{:a;N;/#include ".\\/charmm36-jul2022.ff\\/forcefield.itp"/!ba;s/#include ".\\/charmm36-jul2022.ff\\/forcefield.itp"/&\\n\\n; Include ligand parameters\\n#include "unk.prm"/}' topol.top
sed -i '/\\\[ molecules \\\]/,/Protein_chain_A/ {/Protein_chain_A/ a\\
UNK 1
}' topol.top
```

- _Define the simulation box:_ Position the protein–ligand complex in a simulation box with sufficient volume such that there is at least ~1 nm (or a chosen buffer) between the protein and the box edges. Common practice is to use a **periodic box** (often cubic or truncated octahedral for efficiency). For example, gmx editconf can enlarge the box to the desired dimensions​ (7). An octahedral box is convenient for roughly spherical proteins, as it requires fewer water molecules than a cubic box for the same minimum solute–wall distance.

```
gmx editconf -f complex.gro -o newbox.gro -bt cubic -c -d 1.0
```

- _Solvate the system:_ Fill the box with water molecules (explicit solvent). Using gmx solvate, water molecules (TIP3P model) are added around the protein–ligand, creating a hydration shell (8).

```
gmx solvate -cp newbox.gro -cs spc216.gro -p topol.top -o solv.gro
```

After solvation, the topology is updated to include water molecules.

- _Add ions:_ Replace some water molecules with counter-ions (e.g., Na+ and Cl-) to neutralize the system’s net charge and achieve a realistic salt concentration (e.g., 0.15 M NaCl for physiological conditions). GROMACS’s gmx genion tool can be used after solvation to insert ions at positions of favorable electrostatic potential. This yields a **neutral, solvated, salted system** ready for simulation.

```
cat &lt;<EOF &gt; ions.mdp
; Neighbor searching
cutoff-scheme = Verlet
rlist = 1.1
pbc = xyz
verlet-buffer-tolerance = -1
; Electrostatics
coulombtype = PME
pme-order = 4
fourierspacing = 0.10
rcoulomb = 1.0
; Van der Waals interactions
rvdw = 1.0
; Simulation settings
integrator = steep ; Steepest descent minimization
emtol = 1000.0 ; Convergence criterion (kJ/mol/nm)
nsteps = 50000 ; Maximum number of steps
; Constraints (set to none for minimization)
constraints = none
; Continuation setting
continuation = yes
EOF

gmx grompp -f ions.mdp -c solv.gro -p topol.top -o ions.tpr -maxwarn 1
gmx genion -s ions.tpr -o solv_ions.gro -p topol.top -pname NA -nname CL -neutral -conc 0.15 <<EOF
15
EOF
```

At this point, we have a complete system (protein + ligand + water + ions) described by the CHARMM36 force field.

**Energy minimization:** Prior to MD, perform energy minimization to remove bad contacts or steric clashes introduced during setup. In GROMACS, one typically runs a steepest descent (and sometimes conjugate gradient) minimization (gmx mdrun on a minimization .tpr file) until the maximum force falls below a threshold (e.g., <1000 kJ/mol/nm). This relaxes the system’s high-energy overlaps without imparting any thermal motion​ (9).

```
cat &lt;<EOF &gt; em.mdp
; Energy minimization parameters
integrator = steep ; Algorithm (steepest descent minimization)
emtol = 100.0 ; Stop minimization when the maximum force < 100.0 \[kJ/mol/nm\]
emstep = 0.001 ; Energy step size \[nm\]
nsteps = 10000 ; Maximum number of minimization steps
; Output settings (optional)
nstenergy = 10 ; Write energy to log every 10 steps
; Neighbor searching
cutoff-scheme = Verlet ; Pair list with buffering
nstlist = 20 ; Frequency to update the neighbor list
ns-type = grid ; Make a grid in the box
verlet-buffer-tolerance = 0.005 ; Max allowed error for pair interactions
pbc = xyz ; Periodic boundary conditions
; Electrostatics
coulombtype = PME ; Particle-Mesh-Ewald electrostatics
pme-order = 4 ; Interpolation order for PME
fourierspacing = 0.10 ; Fourier-space grid point spacing
rcoulomb = 1.2 ; Coulomb cut-off distance
; Van der Waals interactions
rvdw = 1.2 ; Lennard-Jones cut-off distance
DispCorr = EnerPres ; Apply dispersion correction for energy & pressure
; Constraints
constraints = none ; No constraints (needed for minimization)
EOF

gmx grompp -f em.mdp -c solv_ions.gro -p topol.top -o em.tpr
gmx mdrun -v -deffnm em
```

**Equilibration:** After minimization, the system is equilibrated in steps:

Often, harmonic position restraints are applied to the heavy atoms of the protein (and sometimes the ligand) during this phase to allow water and ions to relax around the protein without significantly perturbing the protein structure.

```
gmx make_ndx -f unk.gro -o index_unk.ndx << EOF
0 & ! a H\*
q
EOF

gmx genrestr -f unk.gro -n index_unk.ndx -o posre_unk.itp -fc 1000 1000 1000 <<EOF
3
EOF

sed -i '/; Include water topology/{x;s/.\*/\\n; Ligand position restraints\\n#ifdef POSRES\\n#include "posre_unk.itp"\\n#endif\\n/;G}' topol.top
```

**Create index file for the complex**

```
gmx make_ndx -f em.gro -o index.ndx << EOF
1 | 13
q
EOF
```

- _NVT equilibration_ (constant Number, Volume, Temperature): The system is brought to the desired temperature (e.g., 300 K) gradually. Typically, one uses a short MD run (several hundred picoseconds) with the **volume fixed** and a thermostat (e.g., velocity-rescaling or Nose–Hoover) to regulate temperature​ (9).

```
cat &lt;<EOF &gt; nvt.mdp
title = Protein-ligand complex NVT equilibration
define = -DPOSRES ; Position restrain the protein and ligand
; Run parameters
integrator = md ; Leap-frog integrator
nsteps = 50000 ; 2 \* 50000 = 100 ps
dt = 0.002 ; 2 fs timestep
; Output control
nstenergy = 500 ; Save energies every 1.0 ps
nstlog = 500 ; Update log file every 1.0 ps
nstxout-compressed = 1000 ; Save compressed coordinates every 2.0 ps
; Bond parameters
continuation = no ; First dynamics run
constraint_algorithm = lincs ; Holonomic constraints
constraints = h-bonds ; Constrain all hydrogen bonds (all bonds for better stability)
lincs_iter = 1 ; Accuracy of LINCS
lincs_order = 4 ; Also related to accuracy
; Neighbor searching and vdW
cutoff-scheme = Verlet
ns_type = grid ; Search neighboring grid cells
nstlist = 20 ; Update neighbor list every 20 steps
rlist = 1.2 ; Neighbor list cutoff (nm)
vdwtype = cutoff
vdw-modifier = force-switch
rvdw-switch = 1.0
rvdw = 1.2 ; Short-range van der Waals cutoff (nm)
; Electrostatics
coulombtype = PME ; Particle Mesh Ewald for long-range electrostatics
rcoulomb = 1.2 ; Short-range electrostatic cutoff (nm)
pme_order = 4 ; Cubic interpolation order
fourierspacing = 0.16 ; Grid spacing for FFT
; Temperature coupling
tcoupl = V-rescale ; Modified Berendsen thermostat
tc-grps = Protein_UNK Water_and_ions ; Two coupling groups (ensure correct names)
tau_t = 0.1 0.1 ; Time constant (ps)
ref_t = 300 300 ; Reference temperature (K) for each group
; Pressure coupling
pcoupl = no ; No pressure coupling in NVT
; Periodic boundary conditions
pbc = xyz ; 3-D PBC
; Dispersion correction (not used for proteins with the C36 additive FF)
DispCorr = no
; Velocity generation
gen_vel = yes ; Assign velocities from Maxwell distribution
gen_temp = 300 ; Temperature for Maxwell distribution
EOF

gmx grompp -f nvt.mdp -c em.gro -r em.gro -p topol.top -n index.ndx -o nvt.tpr
gmx mdrun -deffnm nvt
```

- _NPT equilibration_ (constant Number, Pressure, Temperature): Next, the simulation is continued with pressure coupling turned on (e.g., Parrinello–Rahman barostat at 1 atm) to adjust the density. Over a few hundred picoseconds, the box volume will fluctuate and stabilize to achieve the target pressure​ (9). Position restraints may be kept initially and gradually released. By the end of NPT equilibration, the system should be at the correct temperature and pressure, with a realistic solvent density, and with the protein–ligand properly solvated.

```
cat &lt;<EOF &gt; npt.mdp
title = Protein-ligand complex NPT equilibration
define = -DPOSRES ; Position restrain the protein and ligand
; Run parameters
integrator = md ; Leap-frog integrator
nsteps = 50000 ; 2 \* 50000 = 100 ps
dt = 0.002 ; 2 fs timestep
; Output control
nstenergy = 500 ; Save energies every 1.0 ps
nstlog = 500 ; Update log file every 1.0 ps
nstxout-compressed = 1000 ; Save compressed coordinates every 2.0 ps
; Bond parameters
continuation = yes ; Continuing from NVT
constraint_algorithm = lincs ; Holonomic constraints
constraints = h-bonds ; Constrain all hydrogen bonds (all bonds for better stability)
lincs_iter = 1 ; Accuracy of LINCS
lincs_order = 4 ; Also related to accuracy
; Neighbor searching and vdW
cutoff-scheme = Verlet
ns_type = grid ; Search neighboring grid cells
nstlist = 20 ; Update neighbor list every 20 steps
rlist = 1.2 ; Neighbor list cutoff (nm)
vdwtype = cutoff
vdw-modifier = force-switch
rvdw-switch = 1.0
rvdw = 1.2 ; Short-range van der Waals cutoff (nm)
; Electrostatics
coulombtype = PME ; Particle Mesh Ewald for long-range electrostatics
rcoulomb = 1.2 ; Short-range electrostatic cutoff (nm)
pme_order = 4 ; Cubic interpolation order
fourierspacing = 0.16 ; Grid spacing for FFT
; Temperature coupling
tcoupl = V-rescale ; Modified Berendsen thermostat
tc-grps = Protein_UNK Water_and_ions ; Two coupling groups (ensure correct names)
tau_t = 0.1 0.1 ; Time constant (ps)
ref_t = 300 300 ; Reference temperature (K) for each group
; Pressure coupling (Improved Settings)
pcoupl = Parrinello-Rahman ; More accurate NPT ensemble
pcoupltype = isotropic ; Uniform scaling of box vectors
tau_p = 2.0 ; Pressure time constant (ps)
ref_p = 1.0 ; Reference pressure (bar)
compressibility = 4.5e-5 ; Isothermal compressibility of water, bar^-1
; Reference coordinate scaling (Consistent with constraints)
refcoord_scaling = all
; Periodic boundary conditions
pbc = xyz ; 3-D PBC
; Dispersion correction (not used for proteins with the C36 additive FF)
DispCorr = no
; Velocity generation
gen_vel = no ; No velocity generation after NVT
EOF

gmx grompp -f npt.mdp -c nvt.gro -t nvt.cpt -r nvt.gro -p topol.top -n index.ndx -o npt.tpr
gmx mdrun -deffnm npt
```

It’s common to monitor properties like temperature, pressure, and density during equilibration (via the GROMACS. edr energy output) to ensure they have stabilized.

**Production run (100 ns MD):** With equilibration complete, we remove any remaining restraints and run the full MD production simulation in the NPT ensemble. In our case, we simulate **100 ns at 300 K and 1 atm** (with periodic boundaries) using a 2 fs time step. Coordinates are saved at regular intervals (e.g., every 10 ps) to produce a trajectory file (.xtc or .trr). During this production phase, the protein and ligand are free to move and interact. Modern GPUs or HPC resources allow 100 ns to be completed in reasonable wall-clock time for a system of this size. (If needed, one could split this into multiple shorter runs with checkpoints, though here we assume one continuous trajectory for simplicity.) By the end, we have a trajectory capturing the detailed motion of the complex over 100 ns.

```
cat &lt;<EOF &gt; md.mdp
title = Protein-ligand complex MD simulation
; Run parameters
integrator = md ; Leap-frog integrator
nsteps = 50000000 ; 2 \* 50000000 \* 0.002 fs = 100 ns
dt = 0.002 ; 2 fs time step
; Output control
nstenergy = 10000 ; Save energies every 20.0 ps
nstlog = 10000 ; Update log file every 20.0 ps
nstxout-compressed = 10000 ; Save coordinates every 20.0 ps
; Bond parameters
continuation = yes ; Continuing from NPT
constraint_algorithm = lincs ; Holonomic constraints
constraints = h-bonds ; Constrain all hydrogeon bonds
lincs_iter = 1 ; Accuracy of LINCS
lincs_order = 4 ; Related to accuracy
; Neighbor searching and vdW
cutoff-scheme = Verlet
ns_type = grid ; Search neighboring grid cells
nstlist = 20 ; Update neighbor list every 20 steps
rlist = 1.2 ; Neighbor list cutoff (nm)
vdwtype = cutoff
vdw-modifier = force-switch
rvdw-switch = 1.0
rvdw = 1.2 ; Short-range van der Waals cutoff (nm)
; Electrostatics
coulombtype = PME ; Particle Mesh Ewald for long-range electrostatics
rcoulomb = 1.2 ; Short-range electrostatic cutoff (nm)
pme_order = 4 ; Cubic interpolation order
fourierspacing = 0.16 ; Grid spacing for FFT
; Temperature coupling
tcoupl = V-rescale ; Modified Berendsen thermostat
tc-grps = Protein_UNK Water_and_ions ; Two coupling groups (ensure correct names)
tau_t = 0.1 0.1 ; Time constant (ps)
ref_t = 300 300 ; Reference temperature (K) for each group
; Pressure coupling
pcoupl = Parrinello-Rahman ; Pressure coupling for NPT
pcoupltype = isotropic ; Uniform scaling of box vectors
tau_p = 2.0 ; Pressure time constant (ps)
ref_p = 1.0 ; Reference pressure (bar)
compressibility = 4.5e-5 ; Isothermal compressibility of water, bar^-1
; Periodic boundary conditions
pbc = xyz ; 3-D PBC
; Dispersion correction (not used for proteins with the C36 additive FF)
DispCorr = no
; Velocity generation
gen_vel = no ; Continuing from NPT equilibration
EOF

gmx grompp -f md.mdp -c npt.gro -t npt.cpt -p topol.top -n index.ndx -o md_0_100.tpr
gmx mdrun -deffnm md_0_100
```

<p align="center">
  <img src="https://github.com/user-attachments/assets/650084ea-34c4-45dc-840e-53c2c82d5910" alt="MD Workflow" width="600"/>
</p>
<p align="center">
  <em>Figure 1:</em> Workflow diagram of a typical MD simulation setup and run (GROMACS). Starting from a PDB structure, a GROMACS topology is generated (<code>pdb2gmx</code>), the complex is placed in a box (<code>editconf</code>), solvated with TIP3P water (<code>solvate</code>), and ions are added. After preparing an input file (<code>grompp</code>), energy minimization and MD runs (<code>mdrun</code>) are conducted, yielding trajectory (.xtc/.trr) and energy (.edr) files for analysis. This pipeline was applied to the protein–ligand system with the CHARMM36 force field.
</p>


**Note on time scale:** 100 ns is sufficient for observing local relaxation (side chain adjustments, loop flexibility) and checking stability of the ligand binding mode. If the ligand is stable, one expects it to remain bound over this timeframe. Large-scale conformational changes or complete unbinding may not occur in 100 ns unless the complex is very unstable. For deeper sampling, multiple 100 ns replicas or longer runs might be used, but for this pipeline we proceed with a single 100 ns trajectory.

**Trajectory Processing Before Analysis**

Before analyzing the trajectory, it is crucial to **post-process the raw trajectory** to avoid artifacts due to periodic boundary conditions (PBC) and molecular diffusion. In an NPT simulation, the protein–ligand complex will freely diffuse in the simulation box and may cross the boundaries of the box. GROMACS uses periodic boundaries, so when part of a molecule crosses one side, it reappears on the opposite side. This can cause the protein or ligand to appear “broken” or jump across the box in the raw trajectory, which would lead to misleading analysis (e.g., an RMSD spike if the protein image changes).

**Periodic boundary correction:** To address this, we generate a **PBC-corrected trajectory** where each molecule is made whole and the complex is kept intact. A typical procedure is:

- **Make molecules whole:** Use gmx trjconv with -pbc whole to ensure each individual molecule (protein, ligand) is continuous, stitching together any broken images.
- **Re-center the complex:** Then use gmx trjconv -pbc mol -center to translate the system so that, for example, the protein (center group) stays in the center of the box, and the ligand remains in the same unit cell as the protein​ (12). The -ur compact flag can be used (for triclinic cells like octahedra) to keep the unit cell compact​ (13). We output a new trajectory (often called _\*\_noPBC.xtc_).

```
gmx trjconv -s md_0_100.tpr -f md_0_100.xtc -o md_pbc.xtc -center -pbc mol -ur compact
# When prompted, select **Protein (group 1)** for both the centering and output groups.
```

After this processing, the protein–ligand complex will no longer diffuse out of the primary cell, and the ligand won’t hop to an adjacent image. All solvent and ions are also properly re-wrapped around the centered complex. We emphasize that **all subsequent analyses will be performed on this PBC-corrected trajectory**​ (14,15). This ensures that calculated distances, RMSDs, etc., reflect real motions rather than artifacts of molecules crossing boundaries.

**Alignment (fitting):** Another preprocessing step is structural alignment of frames. For certain analyses (RMSD, RMSF, PCA), we often remove overall rotation and translation by superimposing each frame onto a reference structure (e.g., the starting structure or average structure). In GROMACS tools like gmx rms and gmx rmsf, this can be done by specifying a fitting group (commonly the protein backbone). In our workflow, when computing RMSD or RMSF, we will fit each frame to the initial minimized structure’s Cα/backbone to focus on internal conformational changes. Fitting prevents global drift from inflating these metrics. (One must **fit after fixing PBC**; performing a fit before making molecules whole can yield incorrect structures. The rule of thumb is: first correct PBC, then fit to reference).

```
gmx trjconv -s md_0_100.tpr -f md_pbc.xtc -o md_0_100_fit.xtc -fit rot+trans
# When prompted, select **Backbone (group 4)** for the centering and **Protein (group 1)** for the output groups.
```

With a **processed, fitted trajectory** in hand, we proceed to analyze various properties. Below is a structured **analysis workflow** covering key analyses and what they tell us about the simulation. We cite GROMACS documentation and literature at each step and suggest visualizations (figures) that one would include in a report or blog post to illustrate the results.

# Analysis Workflow: Key Metrics and Their Interpretation

**Overview of analysis steps:** We will examine structural stability (RMSD), per-residue flexibility (RMSF), overall compactness (radius of gyration), solvent exposure (SASA), and complex conformational landscapes (PCA and free energy landscape). Each analysis provides a different perspective on the MD trajectory:

1. **Root Mean Square Deviation (RMSD) – Structural Stability:** The RMSD measures the average deviation of the structure over time from a reference conformation, typically calculated for the protein’s Cα or backbone atoms. RMSD is a primary indicator of whether the simulation has equilibrated to a stable state. We compute RMSD using GROMACS gmx rms, fitting to the starting structure’s backbone and then calculating the deviation of each frame​ (15). In our pipeline, we separately calculate: (a) **Protein RMSD** (backbone) to monitor protein structural drift, and (b) **Ligand RMSD** (often heavy-atom RMSD of the ligand after fitting the protein) to monitor if the ligand remains in the binding site.

_Interpretation:_ A **stable RMSD plateau** (after an initial rise) indicates the system has settled into an equilibrium ensemble around the reference structure. For example, if the protein backbone RMSD stabilizes around ~0.2–0.3 nm, it suggests the protein’s overall fold remained intact and only small fluctuations occurred. If the ligand’s RMSD (relative to its initial bound pose) stays low (e.g., <0.2 nm), it implies the ligand remains bound in a similar conformation throughout. On the other hand, a continuously rising or fluctuating RMSD might indicate slow conformational drift, insufficient equilibration, or even partial unfolding or dissociation events. It’s common to discard the first few nanoseconds (during which RMSD often climbs as the system relaxes) and consider the RMSD behavior after that as an equilibrium indicator. In summary, **RMSD vs. time plots (Figure 2)** allow us to pinpoint when the system equilibrated and how stable the protein and ligand are in the simulation​ (16).

```
gmx rms -s md_0_100.tpr -f md_0_100_fit.xtc -o rmsd_backbone.xvg -fit rot+trans -tu ns
# When prompted, select **Backbone (group 4)** for both the centering and output groups.

gmx rms -s md_0_100.tpr -f md_0_100_fit.xtc -o rmsd_ligand.xvg -fit rot+trans -tu ns
# When prompted, select **UNK (group 13)** for both the centering and output groups.
```

<p align="center">
  <img src="https://github.com/user-attachments/assets/ac594feb-b923-4831-9c79-089a4af7e7be" alt="MD Workflow" width="600"/>
</p>
<p align="center">
  <em>Figure 2:</em> Example backbone RMSD plot over a 1 ns segment (for illustration). After a brief adjustment (<0.1 ns), the RMSD fluctuates around a steady value (~0.1 nm). In a 100 ns simulation, one would expect an initial rise as the protein–ligand system relaxes, followed by a plateau. A stable RMSD indicates the complex maintains its structural integrity over time.
</p>


1. **Root Mean Square Fluctuation (RMSF) – Per-Residue Flexibility:** RMSF computes the time-averaged fluctuation of each atom (or each residue’s Cα) around its mean position. In practice, we calculate **per-residue RMSF** for the protein to identify which regions are flexible or rigid. Using gmx rmsf (after fitting frames to the reference), we obtain an RMSF value for each residue’s backbone​. This is usually plotted as a curve of RMSF (Å or nm) vs. residue index (or residue name/number). High RMSF peaks indicate residues that move a lot during the simulation (e.g. flexible loops or chain termini), whereas low RMSF indicates stable, rigid regions (e.g. the protein core or secondary structure elements that remained intact).

According to the GROMACS manual, _“RMSF is similar to RMSD, but instead of comparing entire structures, it measures how much each residue fluctuates over the simulation”_​.

Technically, it is the standard deviation of the position of each atom/residue after alignment. Often, results are mapped onto the protein structure (coloring the B-factor field in a PDB by RMSF) to visualize flexible regions.

_Interpretation:_ **RMSF highlights the most dynamic parts of the protein.** For example, one might find that loops at the protein’s surface have RMSF of 0.3–0.5 nm (very flexible), while the well-folded helices or beta-sheets in the core have RMSF <0.1 nm. In a protein–ligand complex, RMSF can identify if the binding site residues are stabilized (lower fluctuation) upon ligand binding, or conversely if some induced fit causes certain regions to become more flexible. If the ligand is stable, often the binding pocket residues show moderate RMSF, whereas if the ligand is unbinding, those residues might spike in flexibility. By comparing RMSF of each residue, one can also pinpoint **potential hinge regions or allosteric sites** that exhibit higher motion. Many publications correlate high RMSF regions with functional motions or unstable segments of the protein​. One could present an **RMSF plot (Figure 3)** and highlight key regions (e.g., “loop 50–55 shows the highest RMSF, indicating a flexible region near the active site”). This per-residue analysis complements RMSD: RMSD tells if the whole structure moved, while RMSF tells _where_ it moved the most.

```
gmx rmsf -s md_0_100.tpr -f md_0_100_fit.xtc -o rmsf_residues.xvg -ox bfactor.pdb
# When prompted, select **Backbone (group 4)** for the output.
```

<p align="center">
  <img src="https://github.com/user-attachments/assets/530d33e1-e565-4630-9021-03e4db12e5e6" alt="MD Workflow" width="600"/>
</p>
<p align="center">
  <em>Figure 3:</em> Example RMSF plot of residue’s atoms. Higher peaks correspond to more flexible regions, often found in loops or chain termini, while lower values indicate stable secondary structure elements. Notable fluctuations around atoms 400–600 and 1600–1800 may reflect flexible loop or surface regions of the protein.
</p>

_Figure 3: Example RMSF plot of residue’s atoms. Higher peaks correspond to more flexible regions, often found in loops or chain termini, while lower values indicate stable secondary structure elements. Notable fluctuations around atoms 400–600 and 1600–1800 may reflect flexible loop or surface regions of the protein._

1. **Radius of Gyration (Rg) – Compactness:** The radius of gyration is the mass-weighted root mean square distance of atoms from their common center of mass. For a protein, Rg reflects how compact or expanded the structure is. We calculate Rg for the protein over time using gmx gyrate (typically on Cα atoms)​. This yields Rg (in nm) as a function of simulation time.

_Interpretation:_ **Rg indicates folding state stability.** If the protein remains in a folded, compact state, Rg will stay relatively constant​. For example, a stable globular protein might have Rg ~1.8 nm ± 0.05 nm throughout the simulation. If the protein **unfolds or expands**, Rg would increase (possibly with more scatter). In our protein–ligand simulation, we expect the protein’s Rg to remain nearly invariant if the protein doesn’t denature; a sudden change might indicate a significant conformational change or partial unfolding. The GROMACS tutorial notes, _“If a protein is stably folded, it will likely maintain a steady Rg; if it unfolds, Rg will change over time”_​. Thus, an **Rg vs. time plot (Figure 4)** serves as a check on overall structural stability, complementing RMSD. In a bound complex, sometimes binding can slightly decrease Rg if the ligand acts like a cross-link that holds the protein tighter, but usually the effect is minor. In our analysis, a roughly constant Rg with small fluctuations implies the protein stays **compact and stable**, which is consistent with a well-behaved simulation​.

```
gmx gyrate -s md_0_100.tpr -f md_0_100_fit.xtc -o gyrate.xvg
# When prompted, select **Protein (group 1)** for the output.
```

<p align="center">
  <img src="https://github.com/user-attachments/assets/1ed6797d-bba8-4c75-975e-c4da7771b23f" alt="MD Workflow" width="600"/>
</p>
<p align="center">
  <em>Figure 4:</em> Example radius of gyration (Rg) over time for a protein. In this short trajectory, Rg fluctuates around ~1.40 nm, indicating a stable, compact structure​. A stable Rg in a 100 ns run would likewise suggest the protein’s tertiary structure remains intact (no unfolding).
</p>

1. **Solvent-Accessible Surface Area (SASA) – Solvent Exposure:** SASA quantifies the surface area of the biomolecule that is accessible to solvent (usually water). It is typically reported in nm² and can be computed with gmx sasa in GROMACS for the whole protein or specific parts (by default it uses a probe radius ~0.14 nm to mimic a water molecule)​ (17,18). For a protein–ligand complex, we are often interested in how the **ligand’s binding affects solvent exposure** of both the ligand and the binding site:
    - **Ligand SASA:** how much of the ligand’s surface is solvent-exposed vs buried in the pocket.
    - **Protein SASA:** how the total protein SASA, or the SASA of the binding site residues, changes upon ligand binding.

We can calculate the **total SASA of the protein** over time, and optionally the SASA of the ligand or certain residue groups (using an index file selection in gmx sasa). The GROMACS tool will output SASA per frame; from there we can plot SASA vs. time or compare average SASA values.

_Interpretation:_ **SASA indicates how “open” or “closed” the structure is.** A _decrease in SASA_ usually means more surface is buried (perhaps due to compaction or ligand binding), while an _increase_ means more surface exposure (perhaps due to partial unfolding or ligand leaving). In our complex, if the ligand stays bound, one expects the **ligand’s SASA to remain low** (the ligand is mostly buried in the binding pocket, with only a small portion exposed to solvent). If the ligand SASA significantly increases at some point, it might indicate the ligand is protruding out or partially dissociating, exposing more of its surface to water. For the protein, binding often causes a slight reduction in SASA because the ligand covers up part of the protein’s surface. For example, in one study the apo protein had SASA ~151.6 nm² versus 149.3 nm² when ligand-bound – this _slight reduction in SASA upon ligand binding indicates the protein’s surface becomes less solvent-exposed, consistent with the ligand filling the pocket and stabilizing that region_​ (19). Monitoring SASA over 100 ns can reveal if the complex remains tightly associated (SASA stays low) or if at times the binding pocket opens up (transient SASA increases). It can also identify **breathing motions** of the protein – e.g., periodic fluctuations in SASA if the protein undergoes expansion/contraction motions. Overall, SASA analysis provides insight into **protein folding stability and ligand burial**. (One must ensure the molecule is whole when computing SASA; if the protein is split across PBC, SASA results would be nonsense​ – another reason we fixed the trajectory first.)


```
gmx sasa -s md_0_100.tpr -f md_0_100_fit.xtc -o sasa.xvg
# When prompted, select **Protein (group 1)** for the output.
```

<p align="center">
  <img src="https://github.com/user-attachments/assets/18d5fc5b-0069-4ad9-972e-40c3728f1cd9" alt="MD Workflow" width="600"/>
</p>
<p align="center">
  <em>Figure 5:</em> Example Solvent Accessible Surface Area (SASA) plot over a 100 ns simulation. SASA fluctuates around ~120 nm², indicating consistent solvent exposure. Stable SASA values typically reflect a folded and compact protein conformation during the trajectory.
</p>


1. **Principal Component Analysis (PCA) – Essential Dynamics:** PCA (also called **essential dynamics** in the MD context) is a statistical technique to identify dominant motion patterns in the trajectory​ (20). In MD, we perform PCA by constructing the covariance matrix of atomic fluctuations (typically for Cα atoms) and diagonalizing it to get eigenvectors (principal components, PCs) and eigenvalues (variance along each PC). GROMACS provides gmx covar to build the covariance matrix and gmx anaeig to analyze the eigenvectors. The first few principal components capture the **largest-amplitude collective motions** of the protein during the simulation. For example, PC1 might describe an opening/closing motion of two domains, PC2 might describe a twisting motion, etc. We usually project the trajectory onto these PCs to see how the system moves along these principal axes (20).

_Approach:_ We fit the trajectory to a reference (to remove overall rotation/translation), then use gmx covar on the protein backbone. After getting eigenvalues, we use gmx anaeig to project the trajectory frames onto, say, the top 2 eigenvectors. This yields principal component coordinates for each frame (PC1 vs time, PC2 vs time), and we can also make a **PC1–PC2 scatter plot** showing the distribution of conformations in the reduced space.

_Interpretation:_ **PCA reveals dominant conformational variants and clustering of states.** A 2D scatter plot of PC1 vs PC2 (Figure 5) often shows clusters of points; each cluster corresponds to a basin of similar conformations. For instance, if the protein sampled two major conformational substates, we might see two clusters in the PCA plot. If the points are all tightly clustered in the center, it means the protein’s motion along those PCs is limited (i.e., very stable structure). If the trajectory gradually drifts in one direction in PC space, it might indicate a slow transition or trend. Sometimes, **the ligand’s binding mode can influence these motions**: e.g., a ligand might restrict a hinge motion, leading to less variance along that mode compared to an apo simulation. PCA is especially useful for **characterizing large-scale motions** (like domain movements, hinge-bending, loop motions) which might not be obvious from RMSD alone​ (21). Typically, one would plot the fraction of variance captured by each PC (the eigenvalue spectrum) to show that the first few PCs capture, say, 80% of the motion. Then, by visualizing extreme projections along PC1/PC2 (e.g., generating representative structures with gmx anaeig -extr), we can describe what those motions physically correspond to. In summary, PCA condenses the trajectory into a few collective motions, allowing an expert to discern if the complex had any noteworthy, concerted movement or if it was essentially just thermal noise. Including a **PCA scatter plot** and mentioning what each cluster/state represents adds depth to the analysis (for an advanced audience, one could further perform clustering or compare PCA of multiple simulations).


```
gmx covar -s md_0_100.tpr -f md_0_100_fit.xtc -o eigenvalues.xvg -v eigenvectors.trr -xpma covar.xpm
gmx anaeig -s md_0_100.tpr -f md_0_100_fit.xtc -v eigenvectors.trr -first 1 -last 2 -proj pc1_pc2.xvg
paste &lt;(awk 'NR&gt;17 {print $1, $2}' pc1_pc2.xvg) &lt;(awk 'NR&gt;17 {print $2}' pc1_pc2.xvg) > PC1PC2.xvg
```

<p align="center">
  <img src="https://github.com/user-attachments/assets/34db793a-f408-4554-bf32-c54e0d280e52" alt="MD Workflow" width="600"/>
</p>
<p align="center">
  <em>Figure 6:</em> Principal Component Analysis (PCA) of protein dynamics. The trajectory is projected onto the first two principal components, revealing major collective motions during the simulation. Clustering or restricted movement in PCA space suggests structural stability, while broad dispersion may indicate conformational transitions.
</p>


1. **Free Energy Landscape (FEL) – Conformational Free Energy Surface:** Building on PCA, we can construct a free energy landscape to map out the thermodynamics of the conformational space. A FEL is often plotted as a 2D heatmap where the axes are two reaction coordinates (for example, PC1 and PC2), and the color represents the free energy (or probability) of the system in that state. High-probability regions correspond to low free energy basins (stable states), and low-probability regions are high-energy (rarely visited). We obtain this by taking the 2D histogram of the trajectory projected onto PC1 and PC2 and applying Boltzmann inversion: ΔG=−kBTln⁡(P)\\Delta G = -k_B T \\ln (P)ΔG=−kB​Tln(P) up to a constant​ (22). GROMACS has the tool gmx sham (stochastic histogram analysis method) that can do this, or one can use Python/R to histogram the PC1–PC2 data and convert to ΔG.

_Interpretation:_ **The FEL highlights metastable states and transition barriers.** In an FEL heatmap (Figure 6), **deep blue wells** (if using a blue-low, red-high color scheme) or **low-energy basins** indicate conformations where the system spent a lot of time – essentially, stable conformational states of the protein–ligand complex. If our 100 ns simulation was stable, we might see essentially one large deep basin, corresponding to the bound state structure (and the surrounding thermal fluctuations). Alternatively, there could be a couple of basins separated by a small ridge – for example, the protein might have two different side-chain packing arrangements or the ligand might bind in two slightly different orientations, each corresponding to a local minimum. The height of the ridge between basins indicates the free energy barrier for transitioning between those states. For a well-behaved simulation in a single bound state, the FEL will show one dominant minimum and maybe some shallow shoulders around it (indicating minor sub-states). In our analysis, we interpret the FEL in conjunction with PCA: **it effectively provides a thermodynamic perspective**, telling us which regions of PCA space are truly distinct states versus just noise. A smooth single-basin FEL would mean the complex has one dominant conformation (good stability), whereas multiple basins might suggest the protein–ligand complex explores different conformations (which could be interesting functionally or could mean insufficient sampling if transitions between them are rare). FELs are an advanced analysis, but for an experienced reader, they convey a lot: stability (depth of the main well), flexibility (breadth of the well), and possible alternate conformations (additional wells). They are especially useful if comparing an apo vs holo simulation to see how a ligand shifts the free energy landscape of the protein.

_Note:_ The absolute values of free energy from a single MD simulation are not meaningful (we set the lowest point to 0), but the relative differences between states are informative. Also, constructing a meaningful FEL requires sufficient sampling along the coordinates of interest – one reason we choose global PCs as the axes, since they capture major motions. If the simulation didn’t sample enough, the FEL might be noisy (which can sometimes be smoothed by longer simulations or combining multiple runs).

```
gmx sham -f PC1PC2.xvg -ls FES.xpm
python xpm2txt.py -f FES.xpm -o fel.dat
```

<p align="center">
  <img src="https://github.com/user-attachments/assets/ded9d321-b39b-46a4-b3ee-cdadd9fa20cb" alt="MD Workflow" width="600"/>
</p>
<p align="center">
  <em>Figure 7:</em> Free Energy Landscape (FEL) of the protein–ligand complex. The surface represents the conformational space projected onto the first two principal components (PC1 and PC2), with free energy depicted along the vertical axis. The deep energy basin indicates a stable structural state, while the red dot marks the global energy minimum sampled during the simulation.
</p>


# Best Practices and Common Pitfalls in MD Analysis

Finally, we highlight some **best practices and pitfalls** to be aware of when analyzing MD trajectories, especially for protein–ligand systems:

- **PBC handling and fitting order:** Always correct for periodic boundary conditions **before** any RMSD/RMSF analysis or distance calculations. If you calculate RMSD on a raw trajectory where the protein jumps out of the box, the values will spike erroneously. The recommended workflow is to make molecules whole and center them (as we did)​), and only then perform any structural alignment (fitting). Fitting a broken molecule can produce incorrect structures. GROMACS tools typically assume a contiguous molecule; as the manual warns, if your molecule is split across PBC, analyses like SASA or Rg will give nonsensical results​. So ensure the input trajectory for analysis has no broken molecules.
- **Choice of reference for RMSD:** The reference structure for RMSD (and initial alignment) should be chosen thoughtfully. Commonly it’s the minimized starting structure or the crystal structure. If the starting structure is high-resolution, use that. In cases where the starting structure is significantly different from the equilibrium (e.g., after a long equilibration), some use the average structure of the trajectory as the reference for RMSD. In our pipeline we used the starting minimized structure. The key is consistency: using the same reference allows comparison across the trajectory. (If comparing multiple simulations, one might use a common reference like the crystal structure for all.) A related tip: if a crystal structure exists, comparing the average simulated structure to the crystal at the end is a good validation – low Cα RMSD to the crystal implies the simulation didn’t drift far from experimental reality.
- **Equilibration time and data production:** Do not trust the early part of the trajectory blindly. It’s standard to **discard an initial equilibration period** (e.g. first 5–10 ns of a 100 ns run) from analyses that assume equilibrium. For instance, if you include the initial rapid RMSD rise in an average, it will skew the RMSD mean higher. Many analyses (RMSF, PCA, FEL) are more meaningful on the equilibrated portion of the trajectory. One can confirm equilibration by observing when properties (RMSD, energy, temperature, etc.) stabilize​. Only after that point should you accumulate statistics for publication-quality results.
- **Sufficient sampling and convergence:** A single 100 ns run may or may not capture all relevant conformational states. If your RMSD is still drifting at 100 ns, or if PCA shows a continuous trajectory trend, it could mean the simulation hasn’t fully converged to sample a stable ensemble. In such cases, extending the simulation or running multiple independent replicates (and then combining analyses) helps ensure that results are robust. For example, PCA especially requires adequate sampling; if you see that the first principal component has a very high variance and the simulation only covered part of that motion, the FEL might suggest more states could be found with longer simulation. Always be cautious in over-interpreting one simulation – check if doubling the simulation time yields similar RMSD/RMSF trends.
- **Interpreting stability vs instability:** A **stable RMSD and energy** does not automatically mean the force field and setup are correct – but it is a minimal requirement. If you observe instability (e.g., the protein unfolding unexpectedly), consider if it’s physical (perhaps the protein is unstable without a partner or membrane) or an artifact (perhaps insufficient equilibration or an issue with the force field/parameters). Ensure that the thermostat and barostat are properly tuned (large fluctuations in temperature or pressure can cause artifacts in structure). Also, check for any spurious restraints or rotation of the simulation box (which can sometimes appear as very smooth periodic oscillations in RMSD – a sign something is off).
- **Analyzing ligand-specific interactions:** In addition to the global metrics we covered, a thorough protein–ligand analysis often includes checking specific interactions: e.g., **hydrogen bond occupancy** between the ligand and key residues (using gmx hbond), **distance monitoring** of important contacts (with gmx distance or custom scripts), and per-residue **interaction energy** decomposition (if energy groups were set, one can extract protein–ligand van der Waals and electrostatic energy components​). These analyses can provide direct insight into what stabilizes the ligand. For example, if a hydrogen bond is present >80% of the time, that’s likely a key interaction. While our pipeline focused on structural and ensemble properties, incorporating these specific analyses is a good practice for a comprehensive study.
- **Visualization and quality control:** It’s often helpful to **visualize the trajectory** (using VMD, PyMOL, or Chimera) alongside quantitative analysis. This can catch odd issues (like a ligand jumping due to PBC imaging which analysis might misinterpret if not corrected). Also, watching the trajectory allows you to correlate plots with actual events (e.g., a jump in RMSD at 50 ns corresponds to the ligand flipping orientation). Make sure to inspect the structure at various time points – for instance, superimpose the start and end structures to see if there’s any notable conformational change.
- **Documentation and reproducibility:** Use GROMACS analysis tools with documented settings (e.g., note which index groups were used for RMSD or RMSF). This ensures that another practitioner could reproduce your analysis. When citing results, as we have done, refer to the GROMACS manual or literature for definitions – for example, clarify that RMSD was backbone-heavy-atom and that the reference was the initial frame, etc. Being explicit about these details prevents confusion, especially for an advanced audience who will scrutinize the methodology.

By following these best practices – preparing a proper simulation, rigorously post-processing, performing diverse analyses, and carefully interpreting the results – one can extract a wealth of information from a 100 ns protein–ligand MD simulation. The CHARMM36 force field and TIP3P water provide a reliable foundation, and GROMACS offers a rich suite of analysis tools to evaluate everything from basic stability to complex conformational landscapes. Armed with these analyses, a researcher can confidently discuss how stable the complex was, which parts of the protein moved, whether the ligand binding pocket stayed intact, and what motions might be functionally relevant – all of which are invaluable for understanding molecular mechanisms or guiding drug design efforts.

# References

1- Yu, T., Sudhakar, N., & Okafor, C. D. (2024). Illuminating ligand-induced dynamics in nuclear receptors through MD simulations. Biochimica et Biophysica Acta (BBA)-Gene Regulatory Mechanisms, 195025.  
2- Huang, J., Rauscher, S., Nawrocki, G., Ran, T., Feig, M., De Groot, B. L., ... & MacKerell Jr, A. D. (2017). CHARMM36m: an improved force field for folded and intrinsically disordered proteins. Nature  methods, 14(1), 71-73.  
3- <https://en.wikipedia.org/wiki/Water_model#:~:text=The%20TIP3P%20model%20implemented%20in,The>  
4- <https://cgenff.com/>  
5- <https://angeloraymondrossi.github.io/workshop/charmm-gromacs-small-organic-molecules-new.html#:~:text=,obtained%20from%20the%20CGenFF%20server>  
6- <https://www.researchgate.net/figure/Flow-diagram-of-MD-analysis-for-studying-the-binding-mechanism-of-inhibitors-epitopes-to_fig4_342129757#:~:text=>  
7- <https://manual.gromacs.org/2025.0/user-guide/flow.html#:~:text=406%20Image%3A%20digraph%20flowchart%20,coordinate%20file%20containing%20molecules%20from>  
8- <https://www.researchgate.net/figure/Flow-diagram-of-MD-analysis-for-studying-the-binding-mechanism-of-inhibitors-epitopes-to_fig4_342129757#:~:text=,In%20the%20next%20step%2C>   
9- <https://www.researchgate.net/figure/Flow-diagram-of-MD-analysis-for-studying-the-binding-mechanism-of-inhibitors-epitopes-to_fig4_342129757#:~:text=simulations%20are%20performed%20in%20an,In%20the%20next%20step%2C>  
10- <https://manual.gromacs.org/2025.0/user-guide/flow.html#:~:text=This%20is%20a%20flow%20chart,402>  
11- <https://manual.gromacs.org/2025.0/user-guide/flow.html#:~:text=URL%3D%22..%2Freference,formats.html%23edr%22%20%5D>  
12- <https://gromacs.bioexcel.eu/t/re-centering-the-protein-in-the-box-after-md-simulation/3853#:~:text=Two%20steps%20%3A>  
13- <https://angeloraymondrossi.github.io/workshop/charmm-gromacs-small-organic-molecules-new.html#:~:text=It%20is%20important%20to%20remember,the%20coordinates%20using%20the%20command>  
14- <https://angeloraymondrossi.github.io/workshop/charmm-gromacs-small-organic-molecules-new.html#:~:text=gmx%20trjconv%20,ur%20compact>  
15- <https://angeloraymondrossi.github.io/workshop/charmm-gromacs-small-organic-molecules-new.html#:~:text=gmx%20trjconv%20,pbc%20mol>  
16- <https://angeloraymondrossi.github.io/workshop/charmm-gromacs-small-organic-molecules-new.html#:~:text=Choose%204%20%28,in%20the%20minimized%2C%20equilibrated%20system>  
17- <https://manual.gromacs.org/current/onlinehelp/gmx-sasa.html#:~:text=,energies%20per%20exposed%20surface%20area>  
18- <https://manual.gromacs.org/current/onlinehelp/gmx-sasa.html#:~:text=The%20average%20and%20standard%20deviation,oa>  
19- Wang, Y., Zhou, Y., & Khan, F. I. (2024). Molecular Insights into Structural Dynamics and Binding Interactions of Selected Inhibitors Targeting SARS-CoV-2 Main Protease. International Journal of Molecular Sciences, 25(24), 13482.  
20- David, Charles C., and Donald J. Jacobs. "Principal component analysis: a method for determining the essential dynamics of proteins." Protein dynamics: Methods and protocols (2014): 193-226.  
21- <https://md-davis.readthedocs.io/en/latest/guides/free_energy_landscapes.html#:~:text=%5C%5B%5CDelta%5Cvarepsilon_i%20%3D%20>  
22- <https://angeloraymondrossi.github.io/workshop/charmm-gromacs-small-organic-molecules-new.html#:~:text=Choose%204%20%28,in%20the%20minimized%2C%20equilibrated%20system>  
