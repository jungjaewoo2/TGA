############################################################################
#
#  Program:         LOESS
#
#  Module:          Makefile
#
#  Purpose:         Top-level Makefile
#
#  Creation date:   March 27, 2006
#
#  Modified:
#
#  Send bug reports, comments or suggestions to kostrun@phys.uconn.edu
#
############################################################################


############################################################################
#
#  > make rlabplus
#  creates the LOESS library for rlabplus and cleans object files afterwards
#
#  > make clean
#  to remove the object files after the library
#
############################################################################

SRC = ./src
G77 = $(F77)
CC = gcc -I../gc
ECHO = echo
RANLIB = ranlib
LOESS = ../libloess.a
AR = ar
ARFLAGS = -r -v

GCINCLUDE = ../../../gc/include/

FSRCFILES = loessf.f

CSRCFILES = loessc.c loess.c predict.c misc.c

OBJFILES = loessc.o loess.o predict.o misc.o loessf.o

rlabplus: clean lib cleanlib

lib:
	@( \
	cd $(SRC); \
	for f in $(FSRCFILES); \
	do \
		$(ECHO) -n "Compiling $$f"; \
		$(G77) $(FFLAGS) -c $$f; \
		$(ECHO)  " Done!"; \
	done ; \
	for c in $(CSRCFILES); \
	do \
		$(ECHO) -n "Compiling $$c"; \
		$(CC) $(CFLAGS) -I$(GCINCLUDE) -c $$c; \
		$(ECHO)  " Done!"; \
	done ; \
	$(AR) $(ARFLAGS) $(LOESS) $(OBJFILES) ;\
	$(RANLIB) $(LOESS) ; \
	cd ..; \
	)

cleanlib:
	(\
	cd $(SRC); rm -f *.o *~ 2>/dev/null; \
	cd .. \
	)

clean: cleanlib
	rm -f *.a *~ 2>/dev/null


