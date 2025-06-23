/* version.h */

/*  This file is a part of RLaB ("Our"-LaB)
    Copyright (C) 1992, 1993, 1994, 1995  Ian R. Searle

    Copyright (C) 2005 Marijan Kostrun, rlabplus project

    This program is free software; you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation; either version 2 of the License, or
    (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program; if not, write to the Free Software
    Foundation, Inc., 675 Mass Ave, Cambridge, MA 02139, USA.

    See the file ./COPYING
 ***********************************************************************/
#define \
  RLAB_MAJOR_VERSION        "3"
#define \
  RLAB_MINOR_VERSION        "1"
#define \
  RLAB_PATCH_MAJOR_VERSION  "1"
#define \
  RLAB_PATCH_MINOR_VERSION  "9"

#define version_string RLAB_MAJOR_VERSION "." RLAB_MINOR_VERSION "." RLAB_PATCH_MAJOR_VERSION "-gcc"

#define RLAB_VERSION_NUM     3

//
// rlab with shared object libraries:
//  libgsl.so, libcblas.so, libblas.so and liblapack.so
//
#define HAVE_RLABPLUS_SO

//
// do we support HDF5 yet?
//
#define HAVE_HDF5_SO

//
// use libcurl for http/https/ftp access within rlab
//
#define HAVE_LIBCURL

//
// use openssl for hash
//
#define HAVE_CRYPTO_HASH

#ifdef __riscos
    static char version_arch[] = "Release 2";
#endif
