############################################################################
#
#  Program:         ODEBIM
#
#  Module:          Makefile
#
#  Purpose:         Top-level Makefile
#
#  Creation date:   March 13, 2006
#
#  Modified:
#
#  Send bug reports, comments or suggestions to kostrun@phys.uconn.edu
#
############################################################################


############################################################################
#
#  > make rlabplus
#  creates the ODEBIM library for rlabplus and cleans object files afterwards
#
#  > make clean
#  to remove the object files after the library
#
############################################################################

SRC = ./src
G77 = $(F77)
ECHO = echo
RANLIB = ranlib
ODEBIM = ../libodebim.a
AR = ar
ARFLAGS = -r -v

#
# hompack: fixpdf.f plus dependencies
#
SRCFILES = bim.f subbim.f bimd.f bimd_sub.f

OBJFILES = bim.o subbim.o bimd.o bimd_sub.o

lib:
	@( \
	cd $(SRC); \
	for f in $(SRCFILES); \
	do \
		$(ECHO) -n "Compiling $$f"; \
		$(G77) $(FFLAGS) -c $$f; \
		$(ECHO)  " Done!"; \
	done ; \
	$(AR) $(ARFLAGS) $(ODEBIM) $(OBJFILES) ;\
	$(RANLIB) $(ODEBIM) ; \
	cd ..; \
	)

cleanlib:
	(\
	cd $(SRC); rm -f *.o *~ 2>/dev/null; \
	cd .. \
	)

clean: cleanlib
	rm -f *.a *~ 2>/dev/null

