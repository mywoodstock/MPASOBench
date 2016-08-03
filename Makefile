# Feel free to compile kernels.F90 without using this Makefile. By default
# (i.e. without -DUSE_ITT_NOTIFY_API) you will not need to compile collect.c.
# The minimum suggested set of flags I would recommend are:
#   -O2 -align array64byte -openmp
# plus: (a) whatever you need to specify the target platform / ISA of interest
# and (b) -DSIMDWIDTH=X where X is the number of complex*16 elements that fit
# into a single SIMD vector (i.e. 4 in the case of AVX512). Widths of 2 and 4
# are supported.

include build/make.inc

CFLAGS = -O3 -g -I.
 
OBJ = kernels.o

###
ifeq ($(VTUNE),1)
FFLAGS += -DUSE_VTUNE
LDFLAGS += $(VTUNE_AMPLIFIER_XE_2016_DIR)/lib64/libittnotify.a
CFLAGS += -I $(VTUNE_AMPLIFIER_XE_2016_DIR)/include
endif

ifeq ($(SDE),1)
FFLAGS += -DUSE_SDE
CFLAGS += -DUSE_SDE
endif


KERNELS = 1

BIN_LIST = $(foreach n, $(KERNELS), $(BIN).$n)
VTUNE_BIN_LIST = $(foreach n, $(KERNELS), $(BIN).vtune.$n)

all : $(BIN_LIST)


vtune : $(VTUNE_BIN_LIST)

$(OBJ) : Makefile

kernels.1.o : kernels.F90 module_itt_sde.o api_itt_sde.o
	$(FC) $(FFLAGS) -D__KERNEL_1 -c $< -o $@ ${CFLAGS}

kernels.2.o : kernels.F90 module_itt_sde.o api_itt_sde.o
	$(FC) $(FFLAGS) -D__KERNEL_2 -c $< -o $@ ${CFLAGS}
                                              
api_itt_sde.o : api_itt_sde.c 
	$(CC) $(CFLAGS) -c $< -o $@

module_itt_sde.o : api_itt_sde.o
	$(FC) $(FFLAGS) -c module_itt_sde.f90 -o $@ ${LDFLAGS}

$(BIN).vtune.% : module_itt_sde.o api_itt_sde.o kernels.%.o
	$(FC) $(FFLAGS) -o $@ $< api_itt_sde.o ${LDFLAGS}

$(BIN).ipm.% : kernels.%.o
	$(FC) $(FFLAGS) -o $@ $< $(IPM)

$(BIN).% : kernels.%.o
	$(FC) $(FFLAGS) -o $@ $< ${LDFLAGS}



.PHONY : clean

clean :
	rm -rf $(BIN_LIST) $(VTUNE_BIN_LIST) *.o *.mod kernels*.optrpt
