/* r_pgplot.h */
/*  This file is a part of RLaB ("Our"-LaB)
   Copyright (C) 1993, 1994  Ian R. Searle

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
#ifndef RLAB_PGPLOT_H
#define RLAB_PGPLOT_H







extern Ent *ent_pgplot_draw_grad_plane (int nargs, Datum args[]);
extern Ent *ent_pgplot_printf (int nargs, Datum args[]);
extern Ent *ent_pgplot_pgqcr (int nargs, Datum args[]);
extern Ent * _pg_cpglen (int nargs, Datum args[]);
extern Ent * _pg_cpgqvp (int nargs, Datum args[]);
extern Ent * _pg_pgqwin (int nargs, Datum args[]);
extern Ent *_pg_color(int nargs, Datum args[]);
extern Ent *_pg_cpgvsiz (int nargs, Datum args[]);
extern Ent *_pg_cpgbeg (int nargs, Datum args[]);
extern Ent *_pg_cpgopen (int nargs, Datum args[]);
extern Ent *_pg_cpgsvp (int nargs, Datum args[]);
extern Ent *_pg_cpgvstd (int nargs, Datum args[]);
extern Ent *_pg_cpgswin (int nargs, Datum args[]);
extern Ent *_pg_cpgwnad (int nargs, Datum args[]);
extern Ent *_pg_cpgslct (int nargs, Datum args[]);
extern Ent *_pg_cpgsubp (int nargs, Datum args[]);
extern Ent *_pg_cpgpanl (int nargs, Datum args[]);
extern Ent *_pg_cpgpap (int nargs, Datum args[]);
extern Ent *_pg_cpgenv (int nargs, Datum args[]);
extern Ent *_pg_cpgline (int nargs, Datum args[]);
extern Ent *_pg_cpgpt (int nargs, Datum args[]);
extern Ent *_pg_cpgimag (int nargs, Datum args[]);
extern Ent *_pg_cpggray (int nargs, Datum args[]);
extern Ent *_pg_cpgwedg (int nargs, Datum args[]);
extern Ent *_pg_cpgconb (int nargs, Datum args[]);
extern Ent *_pg_cpgconl (int nargs, Datum args[]);
extern Ent *_pg_cpgcont (int nargs, Datum args[]);
extern Ent *_pg_cpgconf (int nargs, Datum args[]);
extern Ent *_pg_cpgbox (int nargs, Datum args[]);
extern Ent *_pg_cpgtbox (int nargs, Datum args[]);
extern Ent *_pg_cpgstbg (int nargs, Datum args[]);
extern Ent *_pg_cpgsls (int nargs, Datum args[]);
extern Ent *_pg_cpgslw (int nargs, Datum args[]);
extern Ent *_pg_cpgend (int nargs, Datum args[]);
extern Ent *_pg_cpgclos (int nargs, Datum args[]);
extern Ent *_pg_cpgpage (int nargs, Datum args[]);
extern Ent *_pg_cpgeras (int nargs, Datum args[]);
extern Ent *_pg_cpgetxt (int nargs, Datum args[]);
extern Ent *_pg_cpgtext (int nargs, Datum args[]);
extern Ent *_pg_cpgptxt (int nargs, Datum args[]);
extern Ent *_pg_cpgmtxt (int nargs, Datum args[]);
extern Ent *_pg_cpgscf (int nargs, Datum args[]);
extern Ent *_pg_cpgsch (int nargs, Datum args[]);
extern Ent *_pg_cpgask (int nargs, Datum args[]);
extern Ent *_pg_cpglab (int nargs, Datum args[]);
extern Ent *_pg_cpgupdt (int nargs, Datum args[]);
extern Ent *_pg_cpgsci (int nargs, Datum args[]);
extern Ent *_pg_cpgscir (int nargs, Datum args[]);
extern Ent *_pg_cpgqci (int nargs, Datum args[]);
extern Ent *_pg_cpgqcir (int nargs, Datum args[]);
extern Ent *_pg_cpgbbuf (int nargs, Datum args[]);
extern Ent *_pg_cpgebuf (int nargs, Datum args[]);
extern Ent *_pg_cpgqcol (int nargs, Datum args[]);
extern Ent *_pg_cpgmove (int nargs, Datum args[]);
extern Ent *_pg_cpgdraw (int nargs, Datum args[]);
extern Ent *_pg_cpgsfs (int nargs, Datum args[]);
extern Ent *_pg_cpgshs (int nargs, Datum args[]);
extern Ent *_pg_cpgpoly (int nargs, Datum args[]);
extern Ent *_pg_cpgbin (int nargs, Datum args[]);
extern Ent *_pg_cpghist (int nargs, Datum args[]);
extern Ent *_pg_cpgrnge (int nargs, Datum args[]);
extern Ent *_pg_cpgsave (int nargs, Datum args[]);
extern Ent *_pg_cpgunsa (int nargs, Datum args[]);
extern Ent *_pg_cpgrect (int nargs, Datum args[]);

#endif
