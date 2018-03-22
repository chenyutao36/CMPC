CC=gcc
OPENBLAS_PATH = /opt/OpenBLAS
BLASFEO_PATH = /opt/blasfeo
HPIPM_PATH = /opt/hpipm
QORE_PATH = /home/chen/Documents/Packages/QORE
CFLAGS=-Wall 
OPENBLAS_LIB=$(OPENBLAS_PATH)/lib/libopenblas_haswellp-r0.2.20.a 
HPIPM_LIB =$(HPIPM_PATH)/lib/libhpipm.a
BLASFEO_LIB= $(BLASFEO_PATH)/lib/libblasfeo.a
QORE_LIB = $(QORE_PATH)/bin/libqpsolver_dense.a $(QORE_PATH)/bin/libkktpack_dense.a $(QORE_PATH)/bin/libqpcore.a
QORE_BLAS = $(QORE_PATH)/external/blasfeo/lib

OBJ = main.o

OBJ += common.o
OBJ += casadi_wrapper.o
OBJ += casadi_src.o
OBJ += qp_generation.o
OBJ += qpsolver_hpipm_ocp.o
OBJ += rti_step.o
OBJ += full_condensing.o
OBJ += qpsolver_qore.o

# all: $(OBJ)
all: main

main: $(OBJ)
	$(CC) $(CFLAGS) $(OBJ) -I.. -I$(OPENBLAS_PATH)/include $(OPENBLAS_LIB) -I$(HPIPM_PATH)/include -I$(BLASFEO_PATH)/include $(HPIPM_LIB) -I$(BLASFEO_PATH)/include $(BLASFEO_LIB) -I$(QORE_PATH) -I$(QORE_PATH)/QPSOLVER_DENSE/include $(QORE_LIB) -o main -lm -lpthread -lrt -lblas -L$(QORE_BLAS) -lblasfeo

main.o: main.c
	$(CC) $(CFLAGS) -c main.c

rti_step.o: rti_step.c rti_step.h
	$(CC) $(CFLAGS) -c rti_step.c

qpsolver_hpipm_ocp.o: qpsolver_hpipm_ocp.c qpsolver_hpipm_ocp.h
	$(CC) $(CFLAGS) -I.. -I$(HPIPM_PATH)/include -I$(BLASFEO_PATH)/include -c qpsolver_hpipm_ocp.c

qp_generation.o: qp_generation.c qp_generation.h
	$(CC) $(CFLAGS) -c qp_generation.c

casadi_wrapper.o: casadi_wrapper.c casadi_src.h 
	$(CC) $(CFLAGS) -c casadi_wrapper.c

casadi_src.o: casadi_src.c casadi_src.h 
	$(CC) $(CFLAGS) -c casadi_src.c

common.o: common.c common.h
	$(CC) $(CFLAGS) -c common.c

full_condensing.o: full_condensing.c full_condensing.h
	$(CC) $(CFLAGS) -c full_condensing.c

qpsolver_qore.o: qpsolver_qore.c qpsolver_qore.h
	$(CC) $(CFLAGS) -I.. -I$(QORE_PATH) -I$(QORE_PATH)/QPSOLVER_DENSE/include -c qpsolver_qore.c

clean:
	rm -rf *.o main

