# FDE4NONMEM
The original Fractional Differential Equation Solver subroutine for the implementation of fractional non-linear mixed effects models in NONMEM. Developed by the Pharmacometrics group of the Department of Pharmacy of the National Kapodistrian University of Athens (NKUA). 
The repository consists of the following:
- The FORTRAN subroutine "FDEGL.f90"
- Two example NLME studies of a linear and non-linear fractional model respectively. Each folder contains a NONMEM Control File, a simulated dataset and a .txt file providing the NONMEM command that performs the estimation process for the purpose of demonstrating the use of the fractional solver
- A short manual for the use of the subroutine in the context of a NONMEM study
- The FDE numerical solver using the Grunwald-Letnikov Scheme written in MATLAB
