# GR_code
These are example files for simulating the growth of a cylinder case using a Growth and Remodeling (GR) model. The USERMAT subroutine is written in Fortran, and all input files are designed for use with ANSYS APDL.



To run the input files:



1. First, compile the USERMAT file with a Fortran compiler to export `usermatLib.dll` and `usermatLib.lib`.
2. Then, set the environment variable `ANS_USE_UPF` to the path of `usermatLib.dll` and `usermatLib.lib`.
3. You can then simulate the growth of a cylinder by running the input files `test_user_cylinder.inp` in ANSYS APDL.



Prerequisites:
- Fortran Compiler
- ANSYS APDL version ≥ v221
