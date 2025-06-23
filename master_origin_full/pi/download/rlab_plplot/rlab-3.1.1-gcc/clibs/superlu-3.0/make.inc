############################################################################
#
#  Program:         SuperLU
#
#  Module:          make.inc
#
#  Purpose:         Top-level Definitions
#
#  Creation date:   October 2, 1995
#
#  Modified:	    February 4, 1997  Version 1.0
#		    November 15, 1997 Version 1.1
#		    September 1, 1999 Version 2.0
#
############################################################################
#
# rlab modification by marijan kostrun, IV/2005
# uses optimized blas by Kazushige Goto, installed in /usr/local/lib as
# libblas.so, from http://www.cs.utexas.edu/users/kgoto
# original name of kg's blas: libgoto_p3_256-r0.97.so, or similar
#
###########################################################################

#
#  The name of the libraries to be created/linked to
#
TMGLIB       = tmglib$(PLAT).a
SUPERLULIB   = superlu$(PLAT).a
BLASDEF      = -DUSE_VENDOR_BLAS
BLASLIB      = -lblas -l../gc/gc.a


#
#  The archiver and the flag(s) to use when building archive (library)
#  If your system has no ranlib, set RANLIB = echo.
#
ARCH         = ar
ARCHFLAGS    = cr
RANLIB       = ranlib

CC           = gcc
FORTRAN	     = gfortran
LOADER       = gcc

#
#  C preprocessor defs for compilation (-DNoChange, -DAdd_, or -DUpCase)
#
CDEFS        = -DAdd_
