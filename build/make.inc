# compilers
CC= mpiicc
FC = mpiifort

# flags
FFLAGS = -DKNL -g -convert big_endian -FR -qopenmp #-xMIC-AVX512 #-opt-report=5 #array64byte
LDFLAGS =
LDLIBS = -qopenmp

# additional flags for arch; default SIMDWIDTH; binary name

SIMDWIDTH = 4
FFLAGS += -DSIMDWIDTH=$(SIMDWIDTH)

BIN = kernels.knl
