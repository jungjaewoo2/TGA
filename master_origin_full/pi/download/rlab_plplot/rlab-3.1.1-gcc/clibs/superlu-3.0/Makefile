############################################################################
#
#  Program:         SuperLU
#
#  Module:          Makefile
#
#  Purpose:         Top-level Makefile
#
#  Creation date:   October 2, 1995
#
#  Modified:        February 4, 1997  Version 1.0
#		    November 15, 1997 Version 1.1
#		    September 1, 1999 Version 2.0
#		    October 15, 2003 Version 3.0
#
############################################################################
#
# rlab modification by marijan kostrun, IV/2005
#

include make.inc

all: lib 

lib: superlulib 

clean: cleanlib 

superlulib:
	( cd src; $(MAKE) ; rm *.o 2>/dev/null)


cleanlib:
	rm src/*.o; 
	rm superlu.a

