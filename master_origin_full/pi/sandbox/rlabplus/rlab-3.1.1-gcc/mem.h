/* mem.h */

/*  This file is a part of RLaB ("Our"-LaB)
   Copyright (C) 1995 Ian R. Searle

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
   ********************************************************************** */

#ifndef MEM_H
#define MEM_H

#include "rlab.h"

#ifdef HAVE_GC

    #include "gc.h"

    /* #define GC_DEBUG 1 */

    #define GC_MAIOP GC_malloc_atomic_ignore_off_page

    #ifdef GC_DEBUG
        #define GC_MALLOC GC_malloc_uncollectable
        #define GC_malloc_ignore_off_page GC_MALLOC_ATOMIC
        #define GC_malloc_atomic_ignore_off_page GC_MALLOC_ATOMIC
    #endif  /* GC_DEBUG */

#else

    #define GC_MALLOC malloc
    #define GC_REALLOC realloc
    #define GC_FREE free
    #define GC_malloc_ignore_off_page malloc
    #define GC_malloc_atomic_ignore_off_page malloc
    #define GC_MAIOP malloc

    #define CALLOC calloc

#endif /* HAVE_GC */

#endif /* MEM_H */
