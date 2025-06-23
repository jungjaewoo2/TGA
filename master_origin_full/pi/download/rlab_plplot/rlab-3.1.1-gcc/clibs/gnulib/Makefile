#!/bin/bash
############################################################################
#
#  Send bug reports, comments or suggestions to mkostrun@gmail.com
#
############################################################################

############################################################################
#
#  > make
#
############################################################################

SRC = ./src
GCC = gcc
INCLFLAGS = `pkg-config --cflags glib-2.0` -I./
#INCLFLAGS = -I/usr/include/glib-2.0 -I/usr/lib64/glib-2.0/include

ECHO = echo

GNULIBA = ../gnulib.a
AR = ar
ARFLAGS = -r -v
RANLIB = ranlib

SRCFILES = xmalloc.c gettime.c c-ctype.c  exitfail.c  gnu_parse_datetime.c  parse-datetime.c  timespec.c  xalloc-die.c
OBJFILES = ${SRCFILES:%.c=%.o}

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
	$(AR) $(ARFLAGS) $(GNULIBA) $(OBJFILES) ;\
	$(RANLIB) $(GNULIBA) ; \
	cd ..; \
	)

cleanlib:
	(\
	cd $(SRC); \
	rm *.o 2>/dev/null; \
	cd .. \
	)

clean:
	(\
	rm *.a  *~ 2>/dev/null; \
	cd $(SRC); \
	rm *.o *~ 2>/dev/null;\
	cd .. \
	)


