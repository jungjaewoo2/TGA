############################################################################
#
#  Program:         MINPACK
#
#  Module:          Makefile
#
#  Purpose:         Top-level Makefile
#
#  Creation date:   March 20, 2006
#
#  Modified:
#
#  Send bug reports, comments or suggestions to kostrun@phys.uconn.edu
#
############################################################################


############################################################################
#
#  > make rlabplus
#  creates the MINPACK library for rlabplus and cleans object files afterwards
#
#  > make clean
#  to remove the object files after the library
#
############################################################################

SRC = ./src
G77 = $(F77)
ECHO = echo
RANLIB = ranlib
MINPACK = ../libminpack.a
AR = ar
ARFLAGS = -r -v

#
# minpack: qrsolv.f
#
SRCFILES = qrsolv.f

OBJFILES = qrsolv.o

rlabplus: lib cleanlib

lib:
	@( \
	cd $(SRC); \
	for f in $(SRCFILES); \
	do \
		$(ECHO) -n "Compiling $$f"; \
		$(G77) $(FFLAGS) -c $$f; \
		$(ECHO)  " Done!"; \
	done ; \
	$(AR) $(ARFLAGS) $(MINPACK) $(OBJFILES) ;\
	$(RANLIB) $(MINPACK) ; \
	cd ..; \
	)

cleanlib:
	(\
	cd $(SRC); rm -f *.o *~ 2>/dev/null; \
	cd .. \
	)

clean: cleanlib
	rm -f *.a *~ 2>/dev/null

