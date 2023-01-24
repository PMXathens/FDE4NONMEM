# FDE4NONMEM
This repository contains a Fractional Differential Equation Solver subroutine for the implementation of fractional non-linear mixed effects models in NONMEM. Developed by the Pharmacometrics group of the Department of Pharmacy of the National Kapodistrian University of Athens (NKUA). 
Contributors: Christos Kaikousidis and Aris Dokoumetzidis, (c) 2022.
The repository consists of the following:
- FDEGL.f90: The FORTRAN subroutine
- FDEGL.m: The FDE numerical solver using the Grunwald-Letnikov Scheme written in MATLAB
- A folder containing an example of a Linear FDE model with a Control File and a simulated Dataset
- A folder containing an example of a Non-Linear FDE model with a Control File and a simulated Dataset
- FDEGL_User_Guide.txt: A short user-manual
- An application to a real Diazepam dataset
