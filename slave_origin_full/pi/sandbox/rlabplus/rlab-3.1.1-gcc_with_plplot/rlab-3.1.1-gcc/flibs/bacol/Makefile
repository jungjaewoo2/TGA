############################################################################
#
#  Program:         BACOL
#
#  Module:          Makefile
#
#  Purpose:         Top-level Makefile
#
#  Creation date:   March 28, 2006
#
#  Modified:
#
#  Send bug reports, comments or suggestions to kostrun@phys.uconn.edu
#
############################################################################


############################################################################
#
#  > make rlabplus
#  creates the BACOL library for rlabplus and cleans object files afterwards
#
#  > make clean
#  to remove the object files after the library
#
############################################################################

SRC = ./src
G77 = $(F77)
ECHO = echo
RANLIB = ranlib
BACOL = ../libbacol.a
AR = ar
ARFLAGS = -r -v

#
# bacol: bacol.f and depenencies, most notably dassl et co.
#
SRCFILES = \
bdgefa.f  dgbfa.f   dgbsl.f  bdgesl.f  \
bcrslve.f bcrdcmp.f bcolpnt.f \
bacol.f   ddaslv.f  ddassl.f  ddastp.f  ddajac.f  ddanrm.f \
ddatrp.f  ddaini.f  ddawts.f  \
bbsplvn.f bbsplvd.f eval.f    iniyp.f   res.f     gauleg.f  interv.f  revalu.f \
jac.f     sucstp.f \
caljac.f  divdif.f  idamax.f  meshsq.f  bvalues.f \
calres.f  imtql1.f  pythag.f \
berrest.f imtql2.f  reinit.f  \
errval.f  iniy.f    remesh.f

OBJFILES = \
bdgefa.o  dgbfa.o   dgbsl.o  bdgesl.o  \
bcrslve.o bcrdcmp.o bcolpnt.o \
bacol.o   ddaslv.o  ddassl.o  ddastp.o  ddajac.o  ddanrm.o \
ddatrp.o  ddaini.o  ddawts.o  \
bbsplvn.o bbsplvd.o eval.o    iniyp.o   res.o     gauleg.o  interv.o  revalu.o \
jac.o     sucstp.o \
caljac.o  divdif.o  idamax.o  meshsq.o  bvalues.o \
calres.o  imtql1.o  pythag.o \
berrest.o imtql2.o  reinit.o \
errval.o  iniy.o    remesh.o

rlabplus: clean lib cleanlib

lib:
	@( \
	cd $(SRC); \
	for f in $(SRCFILES); \
	do \
		$(ECHO) -n "Compiling $$f"; \
		$(G77) $(FFLAGS) -c $$f; \
		$(ECHO)  " Done!"; \
	done ; \
	$(AR) $(ARFLAGS) $(BACOL) $(OBJFILES) ;\
	$(RANLIB) $(BACOL) ; \
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

