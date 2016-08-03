# compilers
CC= mpiicc
FC = mpiifort

# flags
FFLAGS = -DKNL -g -O3 -align -qopenmp -xMIC-AVX512 #-opt-report=5 #array64byte
LDFLAGS =
LDLIBS = -qopenmp

# additional flags for arch; default SIMDWIDTH; binary name

SIMDWIDTH = 4
FFLAGS += -DSIMDWIDTH=$(SIMDWIDTH)

BIN = kernels.knl
