############################################################################
#
#  Program:         GENZPAK
#
#  Module:          Makefile
#
#  Purpose:         Top-level Makefile
#
#  Creation date:   March 28, 2006
#
#  Modified:
#
#  Send bug reports, comments or suggestions to kostrun@phys.uconn.edu
#
############################################################################


############################################################################
#
#  > make rlabplus
#  creates the GENZPAK library for rlabplus and cleans object files afterwards
#
#  > make clean
#  to remove the object files after the library
#
############################################################################

SRC = ./src
G77 = $(F77)
ECHO = echo
RANLIB = ranlib
GENZPAK = ../libgenzpak.a
AR = ar
ARFLAGS = -r -v

#
# genzpak: integration in multidimensions
#
SRCFILES = \
dcuhre.f decuhr.f patpack.f simpack.f

OBJFILES = \
dcuhre.o decuhr.o patpack.o simpack.o

rlabplus: clean lib cleanlib

lib:
	@( \
	cd $(SRC); \
	for f in $(SRCFILES); \
	do \
		$(ECHO) -n "Compiling $$f"; \
		$(G77) $(FFLAGS) -c $$f; \
		$(ECHO)  " Done!"; \
	done ; \
	$(AR) $(ARFLAGS) $(GENZPAK) $(OBJFILES) ;\
	$(RANLIB) $(GENZPAK) ; \
	cd ..; \
	)

cleanlib:
	(\
	cd $(SRC); rm -f *.o *~ 2>/dev/null; \
	cd .. \
	)

clean: cleanlib
	rm -f *.a *~ 2>/dev/null

