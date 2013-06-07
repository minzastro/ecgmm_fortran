all:
	f2py --opt=-O3 -c -m ecgmm sorting.f90 lfsr_mod.f90 ecgmm.f90

#This is not used so far
openmp:
	f2py --f90flags="-fopenmp -lgomp" -lgomp --opt=-O3 -c -m ecgmm ecgmm.f90
