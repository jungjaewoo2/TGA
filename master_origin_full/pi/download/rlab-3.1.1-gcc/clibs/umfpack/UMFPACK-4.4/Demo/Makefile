#-------------------------------------------------------------------------------
# compile the UMFPACK demos (for GNU make and original make)
#-------------------------------------------------------------------------------

# UMFPACK Version 4.4, Copyright (c) 2005 by Timothy A. Davis.
# All Rights Reserved.  See ../Doc/License for License.

all: umfpack_di_demo umfpack_zi_demo umfpack_dl_demo umfpack_zl_demo \
	umfpack_simple

include ../Make/Make.include

C = $(CC) $(CFLAGS) $(CONFIG) -I../Include -I../../AMD/Include

INC = ../Include/umfpack.h ../../AMD/Include/amd.h

../Lib/libumfpack.a:
	( cd ../Source ; make )

../../AMD/Lib/libamd.a:
	( cd ../../AMD/Source ; make )

UMFPACK = ../Lib/libumfpack.a ../../AMD/Lib/libamd.a

#-------------------------------------------------------------------------------
# Create the demo programs, run them, and compare the output
#-------------------------------------------------------------------------------

dist: umfpack_di_demo.c umfpack_dl_demo.c umfpack_zi_demo.c umfpack_zl_demo.c

umfpack_simple: umfpack_simple.c $(INC) $(UMFPACK)
	$(C) -o umfpack_simple umfpack_simple.c $(UMFPACK) $(LIB)
	./umfpack_simple

# the GNU rules are simpler:
# umfpack_%_demo.c: umfpack_xx_demo.c umfpack_%_demo.sed
# 	- sed -f umfpack_$*_demo.sed < umfpack_xx_demo.c > umfpack_$*_demo.c
# 
# umfpack_%_demo: umfpack_%_demo.c $(INC) $(UMFPACK)
# 	$(C) -o umfpack_$*_demo umfpack_$*_demo.c $(UMFPACK) $(LIB)
# 	./umfpack_$*_demo > my_umfpack_$*_demo.out

# but do this via brute-force, so we can use just a single Makefile:

# double-precision, int verion:
umfpack_di_demo.c: umfpack_xx_demo.c umfpack_di_demo.sed
	- sed -f umfpack_di_demo.sed < umfpack_xx_demo.c > umfpack_di_demo.c

umfpack_di_demo: umfpack_di_demo.c $(INC) $(UMFPACK)
	$(C) -o umfpack_di_demo umfpack_di_demo.c $(UMFPACK) $(LIB)
	./umfpack_di_demo > my_umfpack_di_demo.out
	- diff umfpack_di_demo.out my_umfpack_di_demo.out

# double-precision, long verion:
umfpack_dl_demo.c: umfpack_xx_demo.c umfpack_dl_demo.sed
	- sed -f umfpack_dl_demo.sed < umfpack_xx_demo.c > umfpack_dl_demo.c

umfpack_dl_demo: umfpack_dl_demo.c $(INC) $(UMFPACK)
	$(C) -o umfpack_dl_demo umfpack_dl_demo.c $(UMFPACK) $(LIB)
	./umfpack_dl_demo > my_umfpack_dl_demo.out
	- diff umfpack_dl_demo.out my_umfpack_dl_demo.out

# complex, int verion:
umfpack_zi_demo.c: umfpack_xx_demo.c umfpack_zi_demo.sed
	- sed -f umfpack_zi_demo.sed < umfpack_xx_demo.c > umfpack_zi_demo.c

umfpack_zi_demo: umfpack_zi_demo.c $(INC) $(UMFPACK)
	$(C) -o umfpack_zi_demo umfpack_zi_demo.c $(UMFPACK) $(LIB)
	./umfpack_zi_demo > my_umfpack_zi_demo.out
	- diff umfpack_zi_demo.out my_umfpack_zi_demo.out

# complex, long verion:
umfpack_zl_demo.c: umfpack_xx_demo.c umfpack_zl_demo.sed
	- sed -f umfpack_zl_demo.sed < umfpack_xx_demo.c > umfpack_zl_demo.c

umfpack_zl_demo: umfpack_zl_demo.c $(INC) $(UMFPACK)
	$(C) -o umfpack_zl_demo umfpack_zl_demo.c $(UMFPACK) $(LIB)
	./umfpack_zl_demo > my_umfpack_zl_demo.out
	- diff umfpack_zl_demo.out my_umfpack_zl_demo.out

#-------------------------------------------------------------------------------
# create a demo program that reads in Harwell/Boeing matrices, and run it
#-------------------------------------------------------------------------------

# the output of "make hb" is in the file umf4.out
hb: $(UMFPACK) umf4 readhb readhb_nozeros readhb_size
	- ./readhb_nozeros < HB/can_24.psa > tmp/A
	- ./readhb_size    < HB/can_24.psa > tmp/Asize
	- ./umf4
	- ./readhb_nozeros < HB/west0067.rua > tmp/A
	- ./readhb_size    < HB/west0067.rua > tmp/Asize
	- ./umf4
	- ./readhb_nozeros < HB/fs_183_6.rua > tmp/A
	- ./readhb_size    < HB/fs_183_6.rua > tmp/Asize
	- ./umf4
	- ./readhb         < HB/fs_183_6.rua > tmp/A
	- ./readhb_size    < HB/fs_183_6.rua > tmp/Asize
	- ./umf4
	- ./readhb         < HB/arc130.rua > tmp/A
	- ./readhb_size    < HB/arc130.rua > tmp/Asize
	- ./umf4
	- ./readhb_nozeros < HB/arc130.rua > tmp/A
	- ./readhb_size    < HB/arc130.rua > tmp/Asize
	- ./umf4
	- ./readhb_nozeros < HB/arc130.rua > tmp/A
	- ./readhb_size    < HB/arc130.rua > tmp/Asize
	- ./umf4 a 1e-6

umf4: umf4.c $(UMFPACK)
	$(C) -o umf4 umf4.c $(UMFPACK) $(LIB)

readhb:  readhb.f
	$(F77) $(F77FLAGS) -o readhb readhb.f $(F77LIB)

readhb_size:  readhb_size.f
	$(F77) $(F77FLAGS) -o readhb_size readhb_size.f $(F77LIB)

readhb_nozeros:  readhb_nozeros.f
	$(F77) $(F77FLAGS) -o readhb_nozeros readhb_nozeros.f $(F77LIB)

#-------------------------------------------------------------------------------
# compile the FORTRAN interface and demo
#-------------------------------------------------------------------------------

fortran: $(UMFPACK) umf4hb.f umf4_f77wrapper.o umf4zhb.f umf4_f77zwrapper.o $(UMFPACK)
	$(F77) $(F77FLAGS) -o umf4hb umf4hb.f umf4_f77wrapper.o \
	    $(UMFPACK) $(LIB)
	- ./umf4hb < HB/west0067.rua > my_umf4hb.out
	- diff my_umf4hb.out umf4hb.out
	$(F77) $(F77FLAGS) -o umf4zhb umf4zhb.f umf4_f77zwrapper.o \
	    $(UMFPACK) $(LIB)
	- ./umf4zhb < HB/qc324.cua > my_umf4zhb.out
	- diff my_umf4zhb.out umf4zhb.out

fortran64: $(UMFPACK) umf4hb64.f umf4_f77wrapper64.o umf4_f77zwrapper64.o $(UMFPACK)
	$(F77) $(F77FLAGS) -o umf4hb64 umf4hb64.f umf4_f77wrapper64.o \
	    $(UMFPACK) $(LIB)
	- ./umf4hb64 < HB/west0067.rua > my_umf4hb64.out
	- diff my_umf4hb64.out umf4hb64.out

umf4_f77wrapper.o: umf4_f77wrapper.c $(INC)
	$(C) -c umf4_f77wrapper.c -o umf4_f77wrapper.o

umf4_f77zwrapper.o: umf4_f77zwrapper.c $(INC)
	$(C) -c umf4_f77zwrapper.c -o umf4_f77zwrapper.o

umf4_f77wrapper64.o: umf4_f77wrapper.c $(INC)
	$(C) -DDLONG -c umf4_f77wrapper.c -o umf4_f77wrapper64.o

umf4_f77zwrapper64.o: umf4_f77zwrapper.c $(INC)
	$(C) -DDLONG -c umf4_f77zwrapper.c -o umf4_f77zwrapper64.o

#-------------------------------------------------------------------------------
# Remove all but the files in the original distribution
#-------------------------------------------------------------------------------

purge: clean
	- $(RM) umfpack_simple a.out
	- $(RM) umfpack_di_demo my_umfpack_di_demo.out
	- $(RM) umfpack_dl_demo my_umfpack_dl_demo.out
	- $(RM) umfpack_zi_demo my_umfpack_zi_demo.out
	- $(RM) umfpack_zl_demo my_umfpack_zl_demo.out
	- $(RM) umf4hb umf4zhb *.umf my_umf4hb.out
	- $(RM) umf4hb64 my_umf4hb64.out my_umf4zhb.out
	- $(RM) umf4 readhb readhb_nozeros readhb_size tmp/*
