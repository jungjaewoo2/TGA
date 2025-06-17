############################################################################
#
#  Program:         CONTIN
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
#  > make lib
#  creates the contin library for rlabplus and cleans object files afterwards
#
#  > make clean
#  to remove the object files after the library
#
############################################################################

SRC = ./src
G77 = $(F77)
ECHO = echo
RANLIB = ranlib
CONTIN = ../libcontin.a
AR = ar
ARFLAGS = -r -v

SRCFILES = banjac.f coqual.f denslv.f denjac.f root.f start.f \
  banslv.f corect.f pitcon.f setstp.f tangnt.f \
  dgedi.f

OBJFILES = banjac.o coqual.o denslv.o denjac.o root.o start.o \
  banslv.o corect.o pitcon.o setstp.o tangnt.o \
  dgedi.o

lib:
	@( \
	cd $(SRC); \
	for f in $(SRCFILES); \
	do \
		$(ECHO) -n "Compiling $$f"; \
		$(G77) $(FFLAGS) -c $$f; \
		$(ECHO)  " Done!"; \
	done ; \
	$(AR) $(ARFLAGS) $(CONTIN) $(OBJFILES) ;\
	$(RANLIB) $(CONTIN) ; \
	cd ..; \
	)

cleanlib:
	(\
	cd $(SRC); rm -f *.o *~ 2>/dev/null; \
	cd .. \
	)

clean: cleanlib
	(\
	rm -f *.a *~ 2>/dev/null; \
	)

