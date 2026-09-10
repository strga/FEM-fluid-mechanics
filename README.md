# FEM Solver for Fluid Mechanics (Potential Flow)

Finite element method solver for 2D irrotational fluid flow, implemented in MATLAB.

## Methods
- Galerkin FEM with weak formulation of the Laplace/Poisson equation
- Isoparametric elements with reference element mapping
- Dirichlet and Neumann boundary conditions
- Verified 2nd order convergence

- ## Cases
- NACA 0012 airfoil — incompressible and compressible potential flow at different angle of attack
- Joukowski airfoil — FEM numerical solution vs. complex-variable analytical solution
- Semicircle domain — convergence tests with Dirichlet, Neumann, Robin BCs
- 1D FEM validation problems

- ## Results
Pressure coefficient (Cp), velocity, and density distributions for each case.
Output figures examples are stored in `FEM_2D/Figures/`.

## Requirements
- MATLAB
- GMSH (for mesh generation — `.geo` and `.msh` files included)
- ParaView for data post-processing

## Run
Open `FEM_2D/solve_potential_flow.m` in MATLAB and run.
For (in)compressible flow (un)comment functions in: `FEM_2D/solve_potential_flow.m` 
For Joukowski airfoil: `Joukowski/Joukowski_profile.m`
