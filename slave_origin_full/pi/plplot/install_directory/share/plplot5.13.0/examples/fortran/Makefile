# examples/fortran/Makefile.examples.in for PLplot
###
### Configure this file with CMake to produce Makefile.examples
###
# Copyright (C) 2002-2016  Alan W. Irwin
# Copyright (C) 2004  Rafael Laboissiere
# Copyright (C) 2004  Andrew Ross
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

SHELL = /bin/bash
Fortran = /usr/bin/gfortran  
EXEEXT = 

PKG_CONFIG_ENV = PKG_CONFIG_PATH="/home/pi/plplot/install_directory/lib/pkgconfig:/home/pi/plplot/install_directory/share/plplot5.13.0/../../lib/pkgconfig/"
FORTRAN_RPATHCMD = -Wl,-rpath -Wl,"/home/pi/plplot/install_directory/lib"

EXECUTABLES_list = \
		 x00f$(EXEEXT) \
		 x01f$(EXEEXT) \
		 x02f$(EXEEXT) \
		 x03f$(EXEEXT) \
		 x04f$(EXEEXT) \
		 x05f$(EXEEXT) \
		 x06f$(EXEEXT) \
		 x07f$(EXEEXT) \
		 x08f$(EXEEXT) \
		 x09f$(EXEEXT) \
		 x10f$(EXEEXT) \
		 x11f$(EXEEXT) \
		 x12f$(EXEEXT) \
		 x13f$(EXEEXT) \
		 x14f$(EXEEXT) \
		 x15f$(EXEEXT) \
		 x16f$(EXEEXT) \
		 x16af$(EXEEXT) \
		 x17f$(EXEEXT) \
		 x18f$(EXEEXT) \
		 x19f$(EXEEXT) \
		 x20f$(EXEEXT) \
		 x21f$(EXEEXT) \
		 x22f$(EXEEXT) \
		 x23f$(EXEEXT) \
		 x24f$(EXEEXT) \
		 x25f$(EXEEXT) \
		 x26f$(EXEEXT) \
		 x27f$(EXEEXT) \
		 x28f$(EXEEXT) \
		 x29f$(EXEEXT) \
		 x30f$(EXEEXT) \
		 x31f$(EXEEXT) \
		 x33f$(EXEEXT)


all: $(EXECUTABLES_list)

clean:
	rm -f $(EXECUTABLES_list)

# target_link_libraries(xfortran plfortrandemolib plplotfortran
# plplotfortranc) Note, the corresponding *.pc file lists both plplotfortran
# and plplotfortranc as PUBLIC libraries.  Note, some Fortran compilers
# (e.g., nagfor) by default do not look for included files (e.g.,
# plfortrandemos.inc included by x21f.f90) in the location of the source
# file when the build is occurring in a different directory so must
# specify an extra -I option for that.
.f90$(EXEEXT):
	$(Fortran)  $< -o $@ $(FORTRAN_RPATHCMD) $(shell $(PKG_CONFIG_ENV) pkg-config  --cflags --libs plplot-fortran) -lplfortrandemolib -I"/home/pi/plplot/install_directory/share/plplot5.13.0"/examples/fortran

.SUFFIXES: .f90 $(EXEEXT)
