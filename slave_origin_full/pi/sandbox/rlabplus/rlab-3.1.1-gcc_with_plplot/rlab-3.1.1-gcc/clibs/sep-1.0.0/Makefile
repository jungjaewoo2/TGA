
OS ?= $(shell sh -c 'uname -s | tr "[A-Z]" "[a-z]"')

SOMAJOR = 0
SOMINOR = 6
SOBUGFIX = 0

ifeq ($(OS), darwin)
SONAME = libsep.dylib
SONAME_MAJOR = libsep.$(SOMAJOR).dylib
SONAME_MAJORMINOR = libsep.$(SOMAJOR).$(SOMINOR).dylib
SONAME_FULL = libsep.$(SOMAJOR).$(SOMINOR).$(SOBUGFIX).dylib
SONAME_FLAG = -install_name
LDPATHENV = DYLD_LIBRARY_PATH
else
ifeq ($(OS), linux)
SONAME = libsep.so
SONAME_MAJOR = libsep.so.$(SOMAJOR)
SONAME_MAJORMINOR = libsep.so.$(SOMAJOR).$(SOMINOR)
SONAME_FULL = libsep.so.$(SOMAJOR).$(SOMINOR).$(SOBUGFIX)
SONAME_FLAG = -soname
LDPATHENV = LD_LIBRARY_PATH
else
$(error OS not yet supported)
endif
endif

INSTALL ?= install
DESTDIR =
PREFIX ?= /usr/local
LIBDIR = $(PREFIX)/lib
INCLUDEDIR = $(PREFIX)/include

CC?=gcc
AR?=ar
CPPFLAGS ?=
LDFLAGS ?=

CPPFLAGS += -Isrc
CFLAGS += -Wall -Wextra -O3  # -Werror
CFLAGS_LIB = $(CFLAGS) -fPIC
LDFLAGS_LIB = $(LDFLAGS) -shared -Wl,$(SONAME_FLAG),$(SONAME_MAJOR)

OBJS = src/analyse.o src/convolve.o src/deblend.o src/extract.o \
       src/lutz.o src/aperture.o src/background.o src/util.o

default: all

src/analyse.o src/convolve.o src/deblend.o src/extract.o src/lutz.o: src/%.o: src/%.c src/extract.h src/sepcore.h src/sep.h
	$(CC) $(CPPFLAGS) $(CFLAGS_LIB) -c src/$*.c -o $@

src/aperture.o: src/aperture.c src/aperture.i src/overlap.h src/sepcore.h src/sep.h
	$(CC) $(CPPFLAGS) $(CFLAGS_LIB) -c src/aperture.c -o $@

src/background.o src/util.o: src/%.o: src/%.c src/sepcore.h src/sep.h
	$(CC) $(CPPFLAGS) $(CFLAGS_LIB) -c src/$*.c -o $@

src/$(SONAME): $(OBJS)
	$(CC) $(LDFLAGS_LIB) $^ -lm -o src/$(SONAME_FULL)
	ln -sf $(SONAME_FULL) src/$(SONAME_MAJORMINOR)
	ln -sf $(SONAME_FULL) src/$(SONAME_MAJOR)
	ln -sf $(SONAME_FULL) src/$(SONAME)

libsep.a: $(OBJS)
	$(AR) rcs libsep.a $^
	rm -f src/*.o 

install: library
	$(INSTALL) -D src/sep.h $(DESTDIR)$(INCLUDEDIR)/sep.h
	$(INSTALL) -D src/$(SONAME_FULL) $(DESTDIR)$(LIBDIR)/$(SONAME_FULL)
	$(INSTALL) -D src/$(SONAME_MAJORMINOR) $(DESTDIR)$(LIBDIR)/$(SONAME_MAJORMINOR)
	$(INSTALL) -D src/$(SONAME_MAJOR) $(DESTDIR)$(LIBDIR)/$(SONAME_MAJOR)
	$(INSTALL) -D src/$(SONAME) $(DESTDIR)$(LIBDIR)/$(SONAME)
	$(INSTALL) -D src/libsep.a $(DESTDIR)$(LIBDIR)/libsep.a

uninstall:
	rm $(DESTDIR)$(INCLUDEDIR)/sep.h
	rm $(DESTDIR)$(LIBDIR)/$(SONAME_FULL)
	rm $(DESTDIR)$(LIBDIR)/$(SONAME_MAJORMINOR)
	rm $(DESTDIR)$(LIBDIR)/$(SONAME_MAJOR)
	rm $(DESTDIR)$(LIBDIR)/$(SONAME)
	rm $(DESTDIR)$(LIBDIR)/libsep.a

test: ctest/test_image
	$(LDPATHENV)=src ctest/test_image data/image.fits data/sep.cat
	ctest/compare.py data/image.cat data/sep.cat

ctest/test_image: ctest/test_image.c src/$(SONAME)
	$(CC) $(CPPFLAGS) $(CFLAGS) $(LDFLAGS) -Lsrc ctest/test_image.c -lm -lsep -o ctest/test_image

clean:
	rm -f src/*.o src/*.a src/libsep.* ctest/test_image

all: src/$(SONAME) libsep.a

.PHONY: all test clean library install uninstall
