These are example files for simulating the growth of a cylinder case using a Growth and Remodeling (GR) model. The USERMAT subroutine is written in Fortran, and all input files are designed for use with ANSYS APDL.

To run the input files:

First, compile the USERMAT file with a Fortran compiler to export usermatLib.dll and usermatLib.lib.
Then, set the environment variable ANS_USE_UPF to the path of usermatLib.dll and usermatLib.lib.
You can then simulate the growth of a cylinder by running the input files test_user_cylinder.inp in ANSYS APDL.
Prerequisites:

Fortran Compiler
ANSYS APDL version ≥ v221

test_user_cylinder.inp: main input file for simulate growth and remodeling for a cylinder case
materail.inp: define each element with a material index to define the fiber direction(axial, circum and normal); define the initial material propreties and state varaiable
mesh.inp: define the mesh size for the cylinder