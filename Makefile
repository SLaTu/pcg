CC = mpicc -O3 -std=c99 -std=gnu11
BSC_CC = $(CC)

MKL_LINK = -L${MKLROOT}/lib/intel64 -lmkl_intel_lp64 -lmkl_core -lmkl_sequential -lpthread
MKL_INC = -I${MKLROOT}/include

LAPACK_LINK = /apps/CBLAS/lib/cblas_LINUX.a -L/apps/LAPACK/3.5.0/GCC/lib -lblas
LAPACK_INC = -I/apps/CBLAS/include

MATH = INTEL_MKL
ifeq ($(MATH), INTEL_MKL)
	MATH_INC = $(MKL_INC)
	MATH_LINK = $(MKL_LINK)
endif

ifeq ($(MATH), LAPACK)
	MATH_INC = $(LAPACK_INC)
	MATH_LINK = $(LAPACK_LINK)
endif

CFLAGS = -c -Wall -L/apps/PAPI/5.6.0/lib/ -lpapi -fopenmp -D$(MATH) $(MATH_INC)
LDFLAGS = -lm -L/apps/PAPI/5.6.0/lib/ -lpapi -fopenmp $(MATH_LINK)

ALL = cg.out

ALLSRC = cg_main.c cg_patterns.c cg_aux.c cg_basics.c

ALLOBJ = $(ALLSRC:.c=.o)

$(ALL) : $(ALLOBJ)
	$(BSC_CC) -o $@ $^ $(LDFLAGS)

%.o : %.c
	$(BSC_CC) $(CFLAGS) -o $@ $<

#For the compatibility
.c.o:
	$(BSC_CC) $(CFLAGS) -o $@ $<

.PHONY: clean

clean:
	rm -f $(ALLOBJ) $(ALL) smpcc_*.c mcc_*.c
