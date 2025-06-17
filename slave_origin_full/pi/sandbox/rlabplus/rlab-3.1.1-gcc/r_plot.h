/* r_plot.h */

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

#ifndef RLAB_PLOT_H
#define RLAB_PLOT_H

extern Ent * ent_plplot_printf(int nargs, Datum args[]);
extern Ent * ent_plplot_plchr(int nargs, Datum args[]);
extern Ent * ent_plplot_plgcol0 (int nargs, Datum args[]);
extern Ent * ent_plplot_plmkstrm (int nargs, Datum args[]);
extern Ent * ent_plplot_pllegend (int nargs, Datum args[]);
extern Ent * ent_plplot_plsetopt (int nargs, Datum args[]);
extern Ent * ent_plplot_plfill (int nargs, Datum args[]);
extern Ent * ent_plplot3_printf (int nargs, Datum args[]);

extern Ent * _plot_plprint (int nargs, Datum args[]);
extern Ent * _plot_plreplot (int nargs, Datum args[]);
extern Ent * _plot_plsstrm (int nargs, Datum args[]);
extern Ent * _plot_plssub (int nargs, Datum args[]);
extern Ent * _plot_plinit (int nargs, Datum args[]);
extern Ent * _plot_plstart (int nargs, Datum args[]);
extern Ent * _plot_plsfile (int nargs, Datum args[]);
extern Ent * _plot_plenv (int nargs, Datum args[]);
extern Ent * _plot_plline (int nargs, Datum args[]);
extern Ent * _plot_plline3 (int nargs, Datum args[]);
extern Ent * _plot_plpoly3 (int nargs, Datum args[]);
extern Ent * _plot_plend (int nargs, Datum args[]);
extern Ent * _plot_plend1 (int nargs, Datum args[]);
extern Ent * _plot_pllab (int nargs, Datum args[]);
extern Ent * _plot_plcol (int nargs, Datum args[]);
extern Ent * _plot_plscolbg (int nargs, Datum args[]);
extern Ent * _plot_plscol0 (int nargs, Datum args[]);
extern Ent * _plot_pllsty (int nargs, Datum args[]);
extern Ent * _plot_plclr (int nargs, Datum args[]);
extern Ent * _plot_plpoin (int nargs, Datum args[]);
extern Ent * _plot_plpoin3 (int nargs, Datum args[]);
extern Ent * _plot_plhist (int nargs, Datum args[]);
extern Ent * _plot_plspage (int nargs, Datum args[]);
extern Ent * _plot_plsstrm (int nargs, Datum args[]);
extern Ent * _plot_pladv (int nargs, Datum args[]);
extern Ent * _plot_plgra (int nargs, Datum args[]);
extern Ent * _plot_pltext (int nargs, Datum args[]);
extern Ent * _plot_plflush (int nargs, Datum args[]);
extern Ent * _plot_plbox (int nargs, Datum args[]);
extern Ent * _plot_plbox3 (int nargs, Datum args[]);
extern Ent * _plot_plvsta (int nargs, Datum args[]);
extern Ent * _plot_plvasp (int nargs, Datum args[]);
extern Ent * _plot_plwind (int nargs, Datum args[]);
extern Ent * _plot_plot3d (int nargs, Datum args[]);
extern Ent * _plot_plmesh (int nargs, Datum args[]);
extern Ent * _plot_plw3d (int nargs, Datum args[]);
extern Ent * _plot_plmtex (int nargs, Datum args[]);
extern Ent * _plot_plspause (int nargs, Datum args[]);
extern Ent * _plot_plwid (int nargs, Datum args[]);
extern Ent * _plot_plptex (int nargs, Datum args[]);
extern Ent * _plot_plfontld (int nargs, Datum args[]);
extern Ent * _plot_plfont (int nargs, Datum args[]);
extern Ent * _plot_plsori (int nargs, Datum args[]);
extern Ent * _plot_plscolor (int nargs, Datum args[]);
extern Ent * _plot_plcont (int nargs, Datum args[]);
extern Ent * _plot_plvpor (int nargs, Datum args[]);
extern Ent * _plot_plvpas (int nargs, Datum args[]);
extern Ent * _plot_plerry (int nargs, Datum args[]);
extern Ent * _plot_plssym (int nargs, Datum args[]);
extern Ent * _plot_plpsty (int nargs, Datum args[]);


#endif /* RLAB_PLOT_H */
