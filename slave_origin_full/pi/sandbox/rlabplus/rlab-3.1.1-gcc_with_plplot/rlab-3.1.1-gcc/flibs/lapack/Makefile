#
#  Top Level Makefile for LAPACK
#  Version 3.0
#  June 30, 1999
#

include make.inc

all: lib

lib: lapacklib

clean: cleanlib
	rm -f $(LAPACKLIB)

lapacklib:
	( cd SRC; $(MAKE) )

cleanlib:
	( cd SRC; $(MAKE) clean )
