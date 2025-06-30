# CPGE_calculation
### This is a code calculating the circular photogalvanic effect (CPGE) tensor based on WannierTools packages.
## Usage guide:
### (1) Put the code file "cpge.f90" into the "src" folder of WannierTools.
### (2) Add the variables: "freqmin", "freqmax", "freqnum", "delta_op" into the namelists "parameters" of the "module para" of the file "module.f90" of WannierTools.
###     freqmin, freqmax: float, minimal and maximal value of optical frequency (eV).
###     freqnum: number of points of frequency.
###     delta_op: The broadening of the Lorentz function (eV).
### (3) Add the following contents into the "readinput.f90" of WannierTools:
        freqmin=freqmin*eV2Hartree
        freqmax=freqmax*eV2Hartree
        delta_op=delta_op*eV2Hartree
### (4) Add "cpge.o" into the variable "OBJ" of the file "Makefile" of WannierTools.
### (5) Recompile the WannierTools.

## If ang questions, please contact me: 2304202694@qq.com
