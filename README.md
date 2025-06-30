# CPGE_calculation
### This is a code calculating the circular photogalvanic effect (CPGE) tensor based on WannierTools packages.
## Usage guide:
### (1) Put the code file "cpge.f90" into the "src" folder of WannierTools.
### (2) Add the variables: "freqmin", "freqmax" (float, minimal and maximal value of optical frequency), "freqnum" (number of points of frequency) into the file "module.f90" of WannierTools.
### (3) Add "cpge.o" into the variable "OBJ" of the file "Makefile" of WannierTools.
### (4) Recompile the WannierTools.
