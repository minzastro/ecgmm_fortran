all:
	f2py --opt=-O3 -c -m ecgmm ecgmm.f90

openmp:
	f2py --f90flags="-fopenmp -lgomp" -lgomp --opt=-O3 -c -m ecgmm ecgmm.f90
