include ../make.inc

#######################################################################
#  This is the makefile to create a library for LAPACK.
#  The files are organized as follows:
#       ALLAUX -- Auxiliary routines called from all precisions
#       SCLAUX -- Auxiliary routines called from both REAL and COMPLEX
#       DZLAUX -- Auxiliary routines called from both DOUBLE PRECISION
#                 and COMPLEX*16
#       SLASRC -- Single precision real LAPACK routines
#       CLASRC -- Single precision complex LAPACK routines
#       DLASRC -- Double precision real LAPACK routines
#       ZLASRC -- Double precision complex LAPACK routines
#
#  The library can be set up to include routines for any combination
#  of the four precisions.  To create or add to the library, enter make
#  followed by one or more of the precisions desired.  Some examples:
#       make single
#       make single complex
#       make single double complex complex16
#  Alternatively, the command
#       make
#  without any arguments creates a library of all four precisions.
#  The library is called
#       lapack.a
#  and is created at the next higher directory level.
#
#  To remove the object files after the library is created, enter
#       make clean
#  On some systems, you can force the source files to be recompiled by
#  entering (for example)
#       make single FRC=FRC
#
#  ***Note***
#  The functions lsame, second, dsecnd, slamch, and dlamch may have
#  to be installed before compiling the library.  Refer to the
#  installation guide, LAPACK Working Note 41, for instructions.
#
#######################################################################
 
ALLAUX = ilaenv.o ieeeck.o lsame.o lsamen.o xerbla.o

DZLAUX = \
	 dlagsy.o dlagge.o dlarot.o dlatms.o dlatmr.o dlatm1.o dlatm2.o dlatm3.o dlaran.o dlarnd.o\
   zlagsy.o zlagge.o zlarot.o zlatms.o zlatmr.o zlatm1.o zlatm2.o zlatm3.o zlarnd.o zlaghe.o


all: double complex16 clean

double: $(DZLAUX) $(ALLAUX)
	$(ARCH) $(ARCHFLAGS) ../$(LAPACKLIB) $(DLASRC) $(ALLAUX) \
	$(DZLAUX)
	$(RANLIB) ../$(LAPACKLIB)

complex16: $(ALLAUX) $(DZLAUX)
	$(ARCH) $(ARCHFLAGS) ../$(LAPACKLIB) $(ZLASRC) $(ALLAUX) \
	$(DZLAUX)
	$(RANLIB) ../$(LAPACKLIB)

$(DZLAUX): $(FRC)

FRC:
	@FRC=$(FRC)

clean:
	rm -f *.o

slamch.o: slamch.f
	$(FORTRAN) $(NOOPT) -c $<
dlamch.o: dlamch.f
	$(FORTRAN) $(FRC)$(NOOPT) -c $<
#	$(FORTRAN) $(OPTS) -c $<

.f.o: 
	$(FORTRAN) $(OPTS) -c $<
