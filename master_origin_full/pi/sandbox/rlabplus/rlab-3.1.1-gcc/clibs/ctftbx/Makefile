#!/bin/bash
# rmake for rlabplus: libtftb.so
# Marijan Kostrun, 2007-2008

# rlab source location:
# starting from 2.2.11 there is an install script which
# also installs headers into /usr/local/include/rlab
############################################################################
#
#  Send bug reports, comments or suggestions to mkostrun@gmail.com
#
############################################################################

############################################################################
#
#  > make
#  creates the TFTBX library for rlabplus and cleans object files afterwards
#
############################################################################

SRC = ./src
GCC = gcc
INCLFLAGS = -I../gc/include -I../rlab/include

ECHO = echo
CTFTBXSO = ../libctftbx.so
SOFLAGS = -static -shared -fPIC

CTFTBXA  = ../libctftbx.a
AR = ar
ARFLAGS = -r -v
RANLIB = ranlib

#
# ctftbx:
#
# SRCFILES = \
#   divers.c af2tfr.c af.c bj.c bud.c cw.c distance.c grd.c hough.c \
#   kernel.c mh.c mhs.c mmce.c page.c pmh.c ppage.c pwv.c ri.c \
#   ridb.c ridbn.c ridh.c ridt.c spwv.c stft.c wv.c zam.c
# OBJFILES = \
#   divers.o af2tfr.o af.o bj.o bud.o cw.o distance.o grd.o hough.o \
#   kernel.o mh.o mhs.o mmce.o page.o pmh.o ppage.o pwv.o ri.o \
#   ridb.o ridbn.o ridh.o ridt.o spwv.o stft.o wv.o zam.o

SRCFILES = tftb.c
OBJFILES = tftb.o

rlabplus: clean obj staticlib cleanlib

obj:
	@( \
	cd $(SRC); \
	for f in $(SRCFILES); \
	do \
    $(ECHO) -n "Compiling $$f"; \
    $(GCC) $(CFLAGS) -c -fPIC $(INCLFLAGS) $$f; \
    $(ECHO)  " Done!"; \
	done ; \
	cd ..; \
  )

staticlib:
	@( \
	cd $(SRC); \
	$(AR) $(ARFLAGS) $(CTFTBXA) $(OBJFILES) ;\
	$(RANLIB) $(CTFTBXA) ; \
	cd ..; \
	)

sharedlib:
	@(\
  cd $(SRC); \
  $(ECHO) -n "Creating shared library"; \
  $(GCC) $(SOFLAGS) -o $(CTFTBXSO) $(OBJFILES); \
  $(ECHO)  " Done!"; \
  cd ..; \
  )

cleanlib:
	(\
  cd $(SRC); rm *.o *~ 2>/dev/null; \
  cd .. \
  )

clean:
	(\
  rm *.so *~ 2>/dev/null; \
  rm *.a  *~ 2>/dev/null; \
  cd $(SRC); rm *.o *~ 2>/dev/null;\
  cd .. \
  )


