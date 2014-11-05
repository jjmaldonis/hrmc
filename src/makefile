


APP = ../rmc

SRC = ran2.f90 globals.f90 read_inputs_v2.f90 model_v2.f90 rmc_functions.f90 mutations.f90 scattering_factors.f90 fft2dc.c slicelib.c islice.c fem1_v2.f90 gr.f90 eam.f90 rmc_v_1_0.f90
OBJ = $(addsuffix .o, $(basename $(SRC)))

USE_LMP = -DUSE_LMP #Uncomment if you want to use lammps
ifdef USE_LMP
LMP_OBJ = /home/jjmaldonis/hrmc/lammps_src/fortran2/LAMMPS.o /home/jjmaldonis/hrmc/lammps_src/fortran2/LAMMPS-wrapper.o
LMP_LIB = /home/jjmaldonis/hrmc/lammps_src/liblammps_odie.a -lmpi_cxx -lstdc++ -lm
LMP_INC = -I/home/jjmaldonis/hrmc/lammps_src/ -I/home/jjmaldonis/hrmc/lammps_src/fortran2
endif

#TIMING = -DTIMING #Uncomment if you want to do timing tests

LIBS=

FC = mpif90
CC = icc


.PHONY: opt debug clean
.SUFFIXES: .f90 .c .o
# Set default ways of creating .o files from .f90 and .c files
%.o: %.f90
	${FC} -c $< ${FCFLAGS} ${LMP_INC}
%.o: %.c
	${CC} -c $< ${CCFLAGS}

# Opt is at the top to make it the default
# Need to define compiler flags for opt / debug
opt: FCFLAGS = -O2 -openmp -fpp $(USE_LMP) $(TIMING)
opt: CCFLAGS = -O2
opt: ${APP}

debug: FCFLAGS = -g -debug all -check all -implicitnone -warn all -openmp -fpp $(USE_LMP)
debug: CCFLAGS = -debug all -g
debug: ${APP}

timing: FCFLAGS = -O2 -openmp -fpp $(USE_LMP) $(TIMING)
timing: CCFLAGS = -O2
timing: TIMING = -DTIMING
timing: ${APP}

${APP}: ${OBJ}
	${FC} ${LMP_OBJ} ${LMP_LIB} ${OBJ} ${FCFLAGS} -o ${APP}

clean:
	rm -f *.mod *.o ${APP}
