############################################################################
#
#  Program:         SLATEC
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
#  creates the SLATEC library for rlabplus and cleans object files afterwards
#
#  > make clean
#  to remove the object files after the library
#
############################################################################

SRC = ./src
G77 = $(F77)
ECHO = echo
RANLIB = ranlib
SLATEC = ../libslatec.a
AR = ar
ARFLAGS = -r -v

#
# slatec: dpolft.f, dpcoef.f, dp1vlu.f
#
SRCFILES = \
dlpdp.f  dp1vlu.f  dwnlit.f  dwnlt2.f  fdump.f   xercnt.f  xerprn.f \
dh12.f   dlsei.f   dpcoef.f  dwnlsm.f  dwnlt3.f  xerhlt.f  xersve.f \
dhfti.f  dlsi.f    dpolft.f  dwnlt1.f  dwnnls.f  i1mach.f  j4save.f \
xermsg.f xgetua.f  davint.f

OBJFILES = \
dlpdp.o  dp1vlu.o  dwnlit.o  dwnlt2.o  fdump.o   xercnt.o  xerprn.o \
dh12.o   dlsei.o   dpcoef.o  dwnlsm.o  dwnlt3.o  xerhlt.o  xersve.o \
dhfti.o  dlsi.o    dpolft.o  dwnlt1.o  dwnnls.o  i1mach.o  j4save.o \
xermsg.o xgetua.o  davint.o

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
	$(AR) $(ARFLAGS) $(SLATEC) $(OBJFILES) ;\
	$(RANLIB) $(SLATEC) ; \
	cd ..; \
	)

cleanlib:
	(\
	cd $(SRC); rm -f *.o *~ 2>/dev/null; \
	cd .. \
	)

clean: cleanlib
	rm -f *.a *~ 2>/dev/null

