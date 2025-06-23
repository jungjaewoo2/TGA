############################################################################
#
#  Program:         HOMPACK
#
#  Module:          Makefile
#
#  Purpose:         Top-level Makefile
#
#  Creation date:   January 19, 2006
#
#  Modified:
#
#  Send bug reports, comments or suggestions to kostrun@phys.uconn.edu
#
############################################################################


############################################################################
#
#  > make rlabplus
#  creates the HOMPACK library for rlabplus and cleans object files afterwards
#
#  > make clean
#  to remove the object files after the library
#
############################################################################

SRC = ./src
G77 = $(F77)
ECHO = echo
RANLIB = ranlib
HOMPACK = ../libhompack.a
AR = ar
ARFLAGS = -r -v

#
# hompack: fixpdf.f plus dependencies
#
SRCFILES = FIXPDF.f FODE.f STEPS.f SINTRP.f ROOT.f\
  FIXPNF.f STEPNF.f ROOTNF.f TANGNF.f\
  FIXPQF.f STEPQF.f ROOTQF.f TANGQF.f UPQRQF.f\
  RHOA.f


OBJFILES = FIXPDF.o  FODE.o STEPS.o SINTRP.o ROOT.o\
  FIXPNF.o STEPNF.o ROOTNF.o TANGNF.o\
  FIXPQF.o STEPQF.o ROOTQF.o TANGQF.o UPQRQF.o\
  RHOA.o

lib:
	@( \
	cd $(SRC); \
	for f in $(SRCFILES); \
	do \
		$(ECHO) -n "Compiling $$f"; \
		$(G77) $(FFLAGS) -c $$f; \
		$(ECHO)  " Done!"; \
	done ; \
	$(AR) $(ARFLAGS) $(HOMPACK) $(OBJFILES) ;\
	$(RANLIB) $(HOMPACK) ; \
	cd ..; \
	)

cleanlib:
	(\
	cd $(SRC); rm -f *.o *~ 2>/dev/null; \
	cd .. \
	)

clean: cleanlib
	rm -f *.a *~ 2>/dev/null

