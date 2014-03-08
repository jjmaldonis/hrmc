


APP = ../rmc
APP_FILE = rmc_v_1_0.f90

SRC = ran2.f90 globals.f90 read_inputs_v2.f90 model_v2.f90 rmc_functions.f90 scattering_factors.f90 fft2dc.c slicelib.c islice.c fem1_v2.f90 gr.f90 eam.f90

OBJ = ran2.o   globals.o   read_inputs_v2.o   model_v2.o   rmc_functions.o   scattering_factors.o   fft2dc.o slicelib.o islice.o fem1_v2.o   gr.o   eam.o

LMP_OBJ = ../lammps_src/fortran2/LAMMPS.o ../lammps_src/fortran2/LAMMPS-wrapper.o
LMP_LIB = ../lammps_src/liblammps_odie.a -lmpi_cxx -lstdc++ -lm
LMP_INC = -I../lammps_src/ -I../lammps_src/fortran2

LIBS=

CFLAGS = -O2 -openmp -fpp
CDEBUGFLAGS = -g -debug all -check all -implicitnone -warn all -O2 -openmp -fpp

FC_OPT = mpif90
CC_OPT = icc

FC_DEBUG = mpif90 -g -debug all -check all -implicitnone -warn all
CC_DEBUG = icc -g  -Wall


.SUFFIXES: .f90 .c .o
%.o: %.mod

all: ${OBJ} ${APP_FILE:.f90=.o} ${APP}

.f90.o:
	${FC_OPT} -c $< ${LMP_INC}

${APP}: ${OBJ} ${APP_FILE:.f90=.o}
	${FC_OPT} ${LMP_OBJ} ${LMP_LIB} ${OBJ} ${APP_FILE:.f90=.o} ${CFLAGS} -o ${APP}

${APP_FILE}: ${OBJ}
	${FC_OPT} -c ${APP_FILE} ${LMP_INC}

.c.o:
	${CC_OPT} -c $<

clean:
	rm -f *.mod *.o ${APP}
