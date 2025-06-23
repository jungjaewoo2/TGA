LIBRARY = libfftpack.a

# F77 =		@F77@
AR = 		ar
RANLIB = 	ranlib
INSTALL =	/usr/bin/install -c
INSTALL_PROGRAM = $(INSTALL)
#FFLAGS =	@FFLAGS@

prefix =        /usr/local
exec_prefix =   ${prefix}
bindir = 	$(exec_prefix)/bin
libdir =        $(exec_prefix)/lib
mand =		$(exec_prefix)/man
mandir = 	$(exec_prefix)/man/man1
buildir =       @BUILD_DIR@/flibs/fftpack

.f.o:
	$(F77) -c $(FFLAGS) $*.f

OBJ = dcfftb.o dcfftf.o dcffti.o dcftb1.o dcftf1.o dcfti1.o dpssb.o dpssb2.o dpssb3.o \
	dpssb4.o dpssb5.o dpssf.o dpssf2.o dpssf3.o dpssf4.o dpssf5.o dradb2.o dradb3.o \
	dradb4.o dradb5.o dradbg.o dradf2.o dradf3.o dradf4.o dradf5.o dradfg.o

$(LIBRARY):	$(OBJ)
		$(AR) r $(LIBRARY) $?
		$(RANLIB) $(LIBRARY)

cleanlib:
		rm -f $(OBJ)

clean:	cleanlib
		rm -f $(LIBRARY)
