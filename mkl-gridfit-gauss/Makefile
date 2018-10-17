OPT=-axAVX -xHost -O2 -B/n/app/openblas/0.2.19/ -B/n/app/openblas/0.2.19/lib

laplace_cpu: main.c cpu.c cpu_band.c Makefile
	icc -vec-report1 $(OPT) cpu_band.c cpu.c main.c -o laplace_cpu_o2 -llapack

#laplace_cpu: main.c cpu.c cpu_band.c Makefile
	#icc -vec-report1 $(OPT) -mkl cpu_band.c cpu.c main.c -o laplace_cpu
