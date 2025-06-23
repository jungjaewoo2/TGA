#-----------------------------------------------------------------------
#            C L A W P A C K  makefile for   RLaB2 Rel.2
#-----------------------------------------------------------------------
#
# use
#    "make lib" to make clawpack.a library for static linking with rlab
#    "make clean" to remove static library and source codes
#
# by Marijan Kostrun, July 14, 2005.

RANLIB = /usr/bin/ranlib

#
# clawpack libraries for now in 1-D
#
DIRS =	1d/lib

lib:
	@for dir in $(DIRS) ;\
	do \
	 (cd $$dir; make libclean; cd ../../) ;\
	done;\
	$(RANLIB) libclaw.a

clean: cleanlib
	rm -f *.a

cleanlib:
	@for dir in $(DIRS) ;\
	do \
	 (cd $$dir; make clean; cd ../../) ;\
	done;\

