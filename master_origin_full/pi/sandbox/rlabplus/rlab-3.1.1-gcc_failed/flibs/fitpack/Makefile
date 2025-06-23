############################################################################
#
#  Program:         FITPAK
#
#  Module:          Makefile
#
#  Purpose:         Top-level Makefile
#
#  Creation date:   June 15, 2007
#
#  Modified:
#
#  Send bug reports, comments or suggestions to kostrun@phys.uconn.edu
#
############################################################################


############################################################################
#
#  > make rlabplus
#  creates the fitpack library for rlabplus and cleans object files afterwards
#
#  > make clean
#  to remove the object files after the library
#
############################################################################

SRC = ./dierckx
G77 = $(F77)
ECHO = echo
RANLIB = ranlib
FITPAK = ../libfit.a
AR = ar
ARFLAGS = -r -v

#
# fitpak
#
SRCFILES = curfit.f clocur.f cocosp.f concur.f cualde.f curev.f curfit.f \
	fourco.f fpader.f fpadno.f fpadpo.f fpback.f fpbacp.f fpbfout.f fpbspl.f \
	fpchec.f fpched.f fpchep.f fpclos.f fpcons.f fpcosp.f fpcsin.f fpcurf.f \
	fpcuro.f fpdeno.f fpdisc.f fpfrno.f fpgivs.f fpinst.f fpintb.f fpknot.f \
	fppara.f fppocu.f fprati.f fprota.f fpseno.f insert.f parcur.f spalde.f \
	spalde.f splder.f splev.f  splint.f sproot.f percur.f fpperi.f \
	concon.f fpcoco.f \
	surfit.f regrid.f fporde.f fprank.f fpsurf.f fpgrre.f fpregr.f \
	parder.f fpbisp.f

OBJFILES = curfit.o clocur.o cocosp.o concur.o cualde.o curev.o curfit.o \
	fourco.o fpader.o fpadno.o fpadpo.o fpback.o fpbacp.o fpbfout.o fpbspl.o \
	fpchec.o fpched.o fpchep.o fpclos.o fpcons.o fpcosp.o fpcsin.o fpcurf.o \
	fpcuro.o fpdeno.o fpdisc.o fpfrno.o fpgivs.o fpinst.o fpintb.o fpknot.o \
	fppara.o fppocu.o fprati.o fprota.o fpseno.o insert.o parcur.o spalde.o \
	spalde.o splder.o splev.o  splint.o sproot.o percur.o fpperi.o \
	concon.o fpcoco.o \
	surfit.o regrid.o fporde.o fprank.o fpsurf.o fpgrre.o fpregr.o \
	parder.o fpbisp.o

lib:
	@( \
	cd $(SRC); \
	for f in $(SRCFILES); \
	do \
		$(ECHO) -n "Compiling $$f"; \
		$(G77) $(FFLAGS) -c $$f; \
		$(ECHO)  " Done!"; \
	done ; \
	$(AR) $(ARFLAGS) $(FITPAK) $(OBJFILES) ;\
	$(RANLIB) $(FITPAK) ; \
	cd ..; \
	)

cleanlib:
	@(\
	cd $(SRC); rm -f *.o *~ 2>/dev/null; \
	cd .. \
	)

clean: cleanlib
	@(\
	rm -f *.a *~ 2>/dev/null; \
	)

