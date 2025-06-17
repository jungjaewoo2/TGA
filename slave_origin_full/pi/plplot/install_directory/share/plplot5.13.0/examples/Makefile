# examples/Makefile.examples.in for PLplot
# configured by our CMake-based build system to examples/Makefile.examples
# in the build tree and installed as examples/Makefile in the install tree.

# Copyright (C) 2004-2017 Alan W. Irwin
#
# This file is part of PLplot.
#
# PLplot is free software; you can redistribute it and/or modify
# it under the terms of the GNU Library General Public License as published by
# the Free Software Foundation; version 2 of the License.
#
# PLplot is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU Library General Public License for more details.
#
# You should have received a copy of the GNU Library General Public License
# along with the file PLplot; if not, write to the Free Software
# Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301  USA

# This Makefile should not be used in the build tree.  It is meant for use
# only in the install tree.  It is first configured, then installed
# (with the name Makefile) in the _installed_ examples directory
# It combines building all the examples that need to be built along
# with testing all examples.

SHELL = /bin/bash
all: c/x01c
c/x01c:
	cd c; $(MAKE)
c_clean:
	cd c; $(MAKE) clean
.PHONY: c_clean
clean: c_clean
x01c.psc: c/x01c
	@echo Generate C results for psc device
	./plplot-test.sh --verbose --front-end=c --device=psc
test_noninteractive: x01c.psc
compare: x01c.psc
# Just for c, c++, and tk directories
test_interactive: c/x01c

#all: ada/x01a
#ada/x01a:
#	cd ada; $(MAKE)
#ada_clean:
#	cd ada; $(MAKE) clean
#.PHONY: ada_clean
#clean: ada_clean
#x01a.psc: ada/x01a
#	@echo Generate Ada results for psc device
#	./plplot-test.sh --verbose --front-end=ada --device=psc
#test_noninteractive: x01a.psc
#compare: x01a.psc

all: c++/x01
c++/x01:
	cd c++; $(MAKE)
cxx_clean:
	cd c++; $(MAKE) clean
.PHONY: cxx_clean
clean: cxx_clean
x01cxx.psc: c++/x01
	@echo Generate C++ results for psc device
	./plplot-test.sh --verbose --front-end=cxx --device=psc
test_noninteractive: x01cxx.psc
compare: x01cxx.psc
# Just for c, c++, and tk directories
test_interactive: c++/x01

#all: d/x01d
#d/x01d:
#	cd d; $(MAKE)
#d_clean:
#	cd d; $(MAKE) clean
#.PHONY: d_clean
#clean: d_clean
#x01d.psc: d/x01d
#	@echo Generate d results for psc device
#	./plplot-test.sh --verbose --front-end=d --device=psc
#test_noninteractive: x01d.psc
#compare: x01d.psc

all: fortran/x01f
fortran/x01f:
	cd fortran; $(MAKE)
fortran_clean:
	cd fortran; $(MAKE) clean
.PHONY: fortran_clean
clean: fortran_clean
x01f.psc: fortran/x01f
	@echo Generate Fortran results for psc device
	./plplot-test.sh --verbose --front-end=fortran --device=psc
test_noninteractive: x01f.psc
compare: x01f.psc

#all: java/x01.class
#java/x01.class:
#	cd java; $(MAKE)
#java_clean:
#	cd java; $(MAKE) clean
#.PHONY: java_clean
#clean: java_clean
#x01j.psc: java/x01.class
#	@echo Generate Java results for psc device
#	env PLPLOT_JAVA_WRAP_DIR="/home/pi/plplot/install_directory/lib/" ./plplot-test.sh --verbose --front-end=java --device=psc
#test_noninteractive: x01j.psc
#compare: x01j.psc

#all: ocaml/x01ocaml
#ocaml/x01ocaml:
#	cd ocaml; $(MAKE)
#ocaml_clean:
#	cd ocaml; $(MAKE) clean
#.PHONY: ocaml_clean
#clean: ocaml_clean
#x01ocaml.psc: ocaml/x01ocaml
#	@echo Generate OCaml results for psc device
#	./plplot-test.sh --verbose --front-end=ocaml --device=psc
#test_noninteractive: x01ocaml.psc
#compare: x01ocaml.psc

#x01o.psc: octave/x01c.m
#	@echo Generate Octave results for psc device
#	./plplot-test.sh --verbose --front-end=octave --device=psc
#test_noninteractive: x01o.psc
#compare: x01o.psc

#x01pdl.psc: perl/x01.pl
#	@echo Generate Perl/PDL results for psc device
#	./plplot-test.sh --verbose --front-end=pdl --device=psc
#test_noninteractive: x01pdl.psc
#compare: x01pdl.psc

#x01p.psc: python/x01
#	@echo Generate Python results for psc device
#	./plplot-test.sh --verbose --front-end=python --device=psc
#test_noninteractive: x01p.psc
#compare: x01p.psc

#x01t.psc: tcl/x01
#	@echo Generate Tcl results for psc device
#	./plplot-test.sh --verbose --front-end=tcl --device=psc
#test_noninteractive: x01t.psc
#compare: x01t.psc

#x01lua.psc: lua/x01.lua
#	@echo Generate Lua results for psc device
#	./plplot-test.sh --verbose --front-end=lua --device=psc
#test_noninteractive: x01lua.psc
#compare: x01lua.psc

#x01c.pdfcairo: c/x01c
#	@echo Generate C results for pdfcairo device
#	./plplot-test.sh --verbose --front-end=c --device=pdfcairo
#test_noninteractive: x01c.pdfcairo

#x01c01.pngcairo: c/x01c
#	@echo Generate C results for pngcairo device
#	./plplot-test.sh --verbose --front-end=c --device=pngcairo
#test_noninteractive: x01c01.pngcairo

#x01c.pscairo: c/x01c
#	@echo Generate C results for pscairo device
#	./plplot-test.sh --verbose --front-end=c --device=pscairo
#test_noninteractive: x01c.pscairo

#x01c01.svgcairo: c/x01c
#	@echo Generate C results for svgcairo device
#	./plplot-test.sh --verbose --front-end=c --device=svgcairo
#test_noninteractive: x01c01.svgcairo

#x01c01.gif : c/x01c
#	@echo Generate C results for gif device
#	./plplot-test.sh --verbose --front-end=c --device=gif
#test_noninteractive: x01c01.gif

#x01c01.jpeg: c/x01c
#	@echo Generate C results for jpeg device
#	./plplot-test.sh --verbose --front-end=c --device=jpeg
#test_noninteractive: x01c01.jpeg

#x01c01.png: c/x01c
#	@echo Generate C results for png device
#	./plplot-test.sh --verbose --front-end=c --device=png
#test_noninteractive: x01c01.png

#x01c.psttfc: c/x01c
#	@echo Generate C results for psttfc device
#	./plplot-test.sh --verbose --front-end=c --device=psttfc
#test_noninteractive: x01c.psttfc

x01c01.svg: c/x01c
	@echo Generate C results for svg device
	./plplot-test.sh --verbose --front-end=c --device=svg
test_noninteractive: x01c01.svg

x01c01.xfig: c/x01c
	@echo Generate C results for xfig device
	./plplot-test.sh --verbose --front-end=c --device=xfig
test_noninteractive: x01c01.xfig

#x01c.pstex: c/x01c
#	@echo Generate C results for pstex device
#	./plplot-test.sh --verbose --front-end=c --device=pstex
#test_noninteractive: x01c.pstex

#x01c.cgm : c/x01c
#	@echo Generate C results for cgm device
#	./plplot-test.sh --verbose --front-end=c --device=cgm
#test_noninteractive: x01c.cgm

#x01c01.bmpqt : c/x01c
#	@echo Generate C results for bmpqt device
#	./plplot-test.sh --verbose --front-end=c --device=bmpqt
#test_noninteractive: x01c01.bmpqt

#x01c01.jpgqt : c/x01c
#	@echo Generate C results for jpgqt device
#	./plplot-test.sh --verbose --front-end=c --device=jpgqt
#test_noninteractive: x01c01.jpgqt

#x01c01.pngqt : c/x01c
#	@echo Generate C results for pngqt device
#	./plplot-test.sh --verbose --front-end=c --device=pngqt
#test_noninteractive: x01c01.pngqt

#x01c01.ppmqt : c/x01c
#	@echo Generate C results for ppmqt device
#	./plplot-test.sh --verbose --front-end=c --device=ppmqt
#test_noninteractive: x01c01.ppmqt

#x01c01.tiffqt : c/x01c
#	@echo Generate C results for tiffqt device
#	./plplot-test.sh --verbose --front-end=c --device=tiffqt
#test_noninteractive: x01c01.tiffqt

#x01c01.svgqt : c/x01c
#	@echo Generate C results for svgqt device
#	./plplot-test.sh --verbose --front-end=c --device=svgqt
#test_noninteractive: x01c01.svgqt

#x01c01.epsqt : c/x01c
#	@echo Generate C results for epsqt device
#	./plplot-test.sh --verbose --front-end=c --device=epsqt
#test_noninteractive: x01c01.epsqt

#x01c01.pdfqt : c/x01c
#	@echo Generate C results for pdfqt device
#	./plplot-test.sh --verbose --front-end=c --device=pdfqt
#test_noninteractive: x01c01.pdfqt

compare:
	@echo Compare non-C language psc and stdout results with corresponding C results
	-./test_diff.sh
	touch compare
test_noninteractive: compare

#all: tk/xtk01$(EXEEXT)
#tk/xtk01$(EXEEXT):
#	cd tk; $(MAKE)
#tk_clean:
#	cd tk; $(MAKE) clean
#.PHONY: tk_clean
#clean: tk_clean
#test_interactive: tk/xtk01$(EXEEXT)

test_interactive:
	./plplot-test-interactive.sh

test_clean:
	rm -f *.psc *.pdfcairo *.pngcairo *.pscairo *.svgcairo \
	*.gif *.jpeg *.png *.psttfc *.svg *.xfig *.pstex* *.*qt *.cgm \
	*_*.txt test.error \
	compare

clean: test_clean

.PHONY : all clean test_noninteractive test_clean test_interactive
