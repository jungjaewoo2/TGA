/* rfile.h */

/*  This file is a part of RLaB ("Our"-LaB)
   Copyright (C) 1994, 1995  Ian R. Searle

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

#ifndef RLAB_RFILE_H
#define RLAB_RFILE_H

extern int rfile_dir (char *dir);

extern void set_search_path _PROTO ((char *value));
extern void set_help_dir _PROTO ((char *value));
extern void set_lib_dir _PROTO ((char *value));
extern void set_pager _PROTO ((char *value));
extern void set_help_pager _PROTO ((char *value));

extern void help _PROTO ((void));
extern void help_name _PROTO ((char *name));
extern void rfile _PROTO ((void));
extern void rfile_load _PROTO ((char *name));
extern void require _PROTO ((char *name));

#endif /* RLAB_RFILE_H */
