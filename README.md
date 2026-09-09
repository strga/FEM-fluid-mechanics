# FEM Solver for Fluid Mechanics (Potential Flow)

Finite element method solver for 2D irrotational fluid flow, implemented in MATLAB.

## Methods
- Galerkin FEM with weak formulation of the Laplace/Poisson equation
- Isoparametric elements with reference element mapping
- Dirichlet and Neumann boundary conditions
- Verified 2nd order convergence

- ## Cases
- NACA 0012 airfoil — incompressible and compressible potential flow at 0°, 4°, 8°
- Joukowski airfoil — FEM numerical solution vs. complex-variable analytical solution
- Semicircle domain — convergence tests with Dirichlet, Neumann, Robin BCs
- 1D FEM validation problems

- ## Results
Pressure coefficient (Cp), velocity, and density distributions for each case.
Output figures examples are stored in `FEM_2D/FEM_2D_figures/`.

## Requirements
- MATLAB
- GMSH (for mesh generation — `.geo` and `.msh` files included)

## Run
Open `FEM_2D/DiscretizePotentialFlow.m` in MATLAB and run.
For compressible flow: `FEM_2D/DiscretizePotentialFlowCompress.m`
For Joukowski airfoil: `Zuk_profil/Zukovskeho_profil.m`
