#
#-----------------------------------------------------------------------
#                     ODRKIT  makefile for RLaB2 Rel. 2
#-----------------------------------------------------------------------
#
# "make" to make the library.
#
# Last Updated:  August 30, 2005.
#

SHELL = /bin/sh
SRC = ./src

G77 = $(F77)
ECHO = echo
RANLIB = ranlib
LIBODR = ../libodr.a
AR = ar
ARFLAGS = -r -v

#
# library objects
#
SRCFILES = d_linpack.f d_odr.f d_odr2.f d_mprec0.f
OBJFILES = d_linpack.o d_odr.o d_odr2.o d_mprec0.o

rlabplus: lib clean

lib:
	@( \
	cd $(SRC); \
	for f in $(SRCFILES); \
	do \
		$(ECHO) -n "Compiling $$f"; \
		$(G77) $(FFLAGS) -c $$f; \
		$(ECHO)  " Done!"; \
	done ; \
	$(AR) $(ARFLAGS) $(LIBODR) $(OBJFILES) ;\
	$(RANLIB) $(LIBODR) ; \
	cd ..; \
	)

cleanlib:
	(cd src; rm -f *~ *.o *.ex .,* fort.* ftn?? *.mat *.pic *.ps *.trace core 2>/dev/null)

clean: cleanlib
	rm -f *.a 2>/dev/null

all: $(OBJ) libodr.a