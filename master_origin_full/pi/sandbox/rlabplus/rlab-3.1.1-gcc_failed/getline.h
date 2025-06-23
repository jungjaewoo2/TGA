/* getline.h: Matrix Dense Real */

/*  This file is a part of RLaB ("Our"-LaB)
   Copyright (C) 2015 M.Kostrun
   Contains reference to functions in getline.c, Copyright (C) 1991 Ian R. Searle, 2015 M. Kostrun

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

#ifndef RLAB_GETLINE_H
#define RLAB_GETLINE_H

//
// getline.c
//
// strings:
extern Ent *Getline (int nargs, Datum args[]);
extern Ent *Strlen (int nargs, Datum args[]);
extern Ent *ent_spinner (int nargs, Datum args[]);
extern Ent *ent_smiley (int nargs, Datum args[]);
extern Ent *ent_string_tolower (int nargs, Datum args[]);
extern Ent *ent_string_toupper (int nargs, Datum args[]);
extern Ent *ent_string_text (int nargs, Datum args[]);
extern Ent *ent_string_blank (int nargs, Datum args[]);
extern Ent *ent_string_ascii (int nargs, Datum args[]);
extern Ent *ent_string_char (int nargs, Datum args[]);
extern Ent *ent_string_substitute (int nargs, Datum args[]);
extern Ent *ent_string_substr (int nargs, Datum args[]);
extern Ent *Scanf   (int nargs, Datum args[]);
extern Ent *Match   (int nargs, Datum args[]);
extern Ent *Grep (int nargs, Datum args[]);
extern Ent *Strsplt (int nargs, Datum args[]);
extern Ent *Findstr (int nargs, Datum args[]);
extern Ent *Strtod (int nargs, Datum args[]);
extern Ent *Lstrip (int nargs, Datum args[]);
extern Ent *Rstrip (int nargs, Datum args[]);
extern Ent *Strip (int nargs, Datum args[]);
extern Ent *ent_string_index (int nargs, Datum args[]);
extern Ent *ent_string_index_pattern (int nargs, Datum args[]);
extern Ent *Argv(int nargs, Datum args[]);
extern Ent *Interactive(int nargs, Datum args[]);
extern Ent *Basename(int nargs, Datum args[]);
extern Ent *ent_Assign(int nargs, Datum args[]);

//
// useful string functions
//
char * process_rlab_string_pattern(char * str);
MDR * mdr_ReadGeneric (FILE * fn, int block_size,
                       int iskiprows, char *delim, int min_line_len, int join_rows, char * join_csp,
                       MDS *comment, MDS *note, MDS *lstrip, MDS *grep,
                       MDR *set_userows, MDR *set_usecols, MDS *start, MDS *stop);


extern Ent *ent_ifelse (int nargs, Datum args[]);

// entropy
extern Ent *ent_iseed (int nargs, Datum args[]);

// terminal
extern Ent *ent_openpty_test (int nargs, Datum args[]);
extern Ent *ent_termcap_clrscr (int nargs, Datum args[]);
extern Ent *ent_termcap_mvcrsr (int nargs, Datum args[]);
extern Ent *ent_termcap_clrpos  (int nargs, Datum args[]);
extern Ent *ent_termcap_colors  (int nargs, Datum args[]);
extern Ent *ent_showkey         (int nargs, Datum args[]);

// system
extern Ent *ent_LsDir (int nargs, Datum args[]);
extern Ent *ent_Stat (int nargs, Datum args[]);
extern Ent *ent_Getpid (int nargs, Datum args[]);
extern Ent *ent_Fnctl (int nargs, Datum args[]);

// time functions
extern Ent *ent_Clock (int nargs, Datum args[]);
extern Ent *ent_time (int nargs, Datum args[]);
extern Ent *ent_gmtime (int nargs, Datum args[]);
extern Ent *ent_Seconds (int nargs, Datum args[]);
extern Ent *ent_strftime (int nargs, Datum args[]);
extern Ent *ent_Etime (int nargs, Datum args[]);
extern Ent *ent_dayofweek (int nargs, Datum args[]);
extern Ent *ent_dayofyear (int nargs, Datum args[]);
extern Ent *ent_dstr2time (int nargs, Datum args[]);




























#endif /* RLAB_GETLINE_H */
