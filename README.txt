Finite Element Analysis and Davidson Eigensolver for Fluid-Solid Phononic Crystals

This repository contains the MATLAB code accompanying the paper:
"An Efficient Structure-Preserving Davidson Algorithm for Fluid-Solid Phononic Crystals".

Overview
Computes eigenvalues, angular frequencies, and modes of three-dimensional periodic
fluid-solid phononic crystals at a specified Bloch wave vector, and records
convergence history, residuals, and computation times. Supports the SC, AlWater,
and SC_Quad models, with initialization using random vectors, prolongated coarse-grid
subsystem modes, fine-grid subsystem modes, or coarse-grid coupled modes.

Usage
In MATLAB, set the current folder to the project root and run:
    SC_Quad_sweep_width_and_all_InitType
Select the model using the model variable in the entry script. Set the wave-vector
index, number of eigenvalues, initialization method, solver tolerances, and other
parameters in separated_Davidson/Initialize_parameters.m.
Each run computes one wave-vector point for one model. For parameter sweeps,
adjust the settings and repeat the run.
Results are saved to data_paper/result_<model>.mat. The result structure contains
eigenvalues, angular frequencies, convergence history, initialization details,
and performance statistics.

Main Functions
- Physical_settings / Physical_settings_AlWater: Define material properties,
  the lattice, and the wave-vector path.
- readGmsh: Read mesh data.
- Assembly_mass_and_stiffness_matrix: Assemble solid and fluid mass and stiffness
  matrices, together with the coupling matrix.
- Generate_Bloch_T_Matrix_fast: Construct reduction matrices for Bloch periodic
  boundary conditions.
- Initialize_parameters / davidson_opts_from_Opt: Configure computation parameters
  and convert them into solver options.
- use_FSI_prolongation_init_time / use_FSI_prolongation_init_time_coupled:
  Construct initial subspaces from subsystem or coupled modes, using Lanczos
  and build_FSI_prolongation_init as needed.
- solve_FSI_separated_Davidson_time_lock: Main separated-subspace Davidson solver
  with mode locking.
- solve_reduced_qep / evaluate_candidates / solve_shifted_system_time:
  Solve the projected quadratic eigenvalue problem, evaluate candidate modes,
  and solve correction equations.
- Prec_setup_solid_ict / Prec_setup_fluid_ict / Prec_setup_Ks_ict:
  Construct preconditioners.
- build_metrics: Summarize solver and initialization performance metrics.
- main_eigs: Function file used for comparison with MATLAB's built-in eigs function.
