#
# Makefile for 1d clawpack code
#
LFLAGS =  
F77 = gfortran -c
AR	= /usr/bin/ar -rcv
RANLIB	= /usr/bin/ranlib

.f.o: ; $(F77) $(FFLAGS) $*.f -o $*.o

LIBOBJECTS = \
  bc1.o \
  setprob.o \
  setaux.o \
  claw1.o \
  step1.o \
  copyq1.o \
  inlinelimiter.o \
  limiter.o \
  philim.o \
  out1.o

SOURCES = \
  bc1.f \
  setprob.f \
  setaux.f \
  claw1.f \
  step1.f \
  copyq1.f \
  inlinelimiter.f \
  limiter.f \
  philim.f \
  out1.f

lib: $(LIBOBJECTS)
	$(AR) ../../libclaw.a $(LIBOBJECTS)

clean:
	rm -f *.o 

libclean: lib clean

### DO NOT remove this line - make depends on it ###
