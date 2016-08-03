# Feel free to compile kernels.F90 without using this Makefile. By default
# (i.e. without -DUSE_ITT_NOTIFY_API) you will not need to compile collect.c.
# The minimum suggested set of flags I would recommend are:
#   -O2 -align array64byte -openmp
# plus: (a) whatever you need to specify the target platform / ISA of interest
# and (b) -DSIMDWIDTH=X where X is the number of complex*16 elements that fit
# into a single SIMD vector (i.e. 4 in the case of AVX512). Widths of 2 and 4
# are supported.

#include Make.inc.intel
include make.inc

#VTUNE_AMPLIFIER_XE_2016_DIR=/opt/intel/vtune/2016u3.external.vtknl/vtune_amplifier_xe

CFLAGS = -g -I $(VTUNE_AMPLIFIER_XE_2016_DIR)/include
 
OBJ = kernels.o

###
ifeq ($(VTUNE),1)
FFLAGS += -DUSE_VTUNE
LDFLAGS += -L$(VTUNE_AMPLIFIER_XE_2016_DIR)/lib64/libittnotify.a
endif

ifeq ($(SDE),1)
FFLAGS += -DUSE_SDE
CFLAGS += -DUSE_SDE
endif


KERNELS = 1 2

BIN_LIST = $(foreach n, $(KERNELS), $(BIN).$n)
VTUNE_BIN_LIST = $(foreach n, $(KERNELS), $(BIN).vtune.$n)

all : $(BIN_LIST)


vtune : $(VTUNE_BIN_LIST)

$(OBJ) : Makefile

kernels.1.o : kernels.F90
	$(FC) $(FFLAGS) -D__KERNEL_1 -c $< -o $@ ${LDFLAGS}

kernels.2.o : kernels.F90
	$(FC) $(FFLAGS) -D__KERNEL_2 -c $< -o $@  ${LDFLAGS}
                                              
api_itt_sde.o : api_itt_sde.c 
	$(CC) $(CFLAGS) -c $< -o $@

module_itt_sde.o : api_itt_sde.o
	$(FC) $(FFLAGS) -c module_itt_sde.f90 -o $@

$(BIN).vtune.% : module_itt_sde.o kernels.%.o
	$(FC) $(FFLAGS) -o $@ $? api_itt_sde.o $(VTUNE_AMPLIFIER_XE_2016_DIR)/lib64/libittnotify.a

$(BIN).ipm.% : kernels.%.o
	$(FC) $(FFLAGS) -o $@ $< $(IPM)

$(BIN).% : kernels.%.o
	$(FC) $(FFLAGS) -o $@ $< 



.PHONY : clean

clean :
	rm -rf $(BIN_LIST) $(VTUNE_BIN_LIST) *.o *.mod kernels*.optrpt
