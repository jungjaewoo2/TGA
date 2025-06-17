/* rimagemagick.h */

/* This file is a part of rlabplus
   Copyright (C) 2013  M. Kostrun

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

extern Ent *ent_im_wand    (int nargs, Datum args[]);
extern Ent *ent_im_init    (int nargs, Datum args[]);
extern Ent *ent_im_clear   (int nargs, Datum args[]);
extern Ent *ent_im_close   (int nargs, Datum args[]);
extern Ent *ent_im_clone   (int nargs, Datum args[]);
extern Ent *ent_im_property   (int nargs, Datum args[]);
extern Ent *ent_im_options    (int nargs, Datum args[]);
extern Ent *ent_im_artifacts  (int nargs, Datum args[]);
extern Ent *ent_im_read_image_to_wand    (int nargs, Datum args[]);
extern Ent *ent_im_wand_write_image (int nargs, Datum args[]);
extern Ent *ent_im_display_wand    (int nargs, Datum args[]);
extern Ent *ent_imlib_display_image(int nargs, Datum args[]);
extern Ent *ent_im_wand_iterate_images(int nargs, Datum args[]);
extern Ent *ent_im_exit(int nargs, Datum args[]);
extern Ent *ent_im_wand_apply(int nargs, Datum args[]);
extern Ent *ent_im_wand_join(int nargs, Datum args[]);
extern Ent *ent_im_wand_append(int nargs, Datum args[]);
extern Ent *ent_im_wand_combine(int nargs, Datum args[]);
extern Ent *ent_im_wand_diff(int nargs, Datum args[]);
extern Ent *ent_im_wand_deconstruct(int nargs, Datum args[]);
extern Ent *ent_im_wand_distort(int nargs, Datum args[]);

// extern Ent *ent_im_wand_func(int nargs, Datum args[]);
// extern Ent *ent_im_wand_eval(int nargs, Datum args[]);
// extern Ent *ent_im_wand_expr(int nargs, Datum args[]);
extern Ent *ent_im_wand_f(int nargs, Datum args[]);

extern Ent *ent_im_wand_clut(int nargs, Datum args[]);
// extern Ent *ent_im_wand_image_constitute(int nargs, Datum args[]);
extern Ent *ent_im_wand_image_pixels(int nargs, Datum args[]);
extern Ent *ent_im_wand_image_specs(int nargs, Datum args[]);

Bltin rlab_im_bltin[] = {
  {BLTIN, "wand",     ent_im_wand},
  {BLTIN, "init",     ent_im_init},
  {BLTIN, "clear",    ent_im_clear},
  {BLTIN, "close",    ent_im_close},
  {BLTIN, "clone",    ent_im_clone},
  {BLTIN, "prop",     ent_im_property},
  {BLTIN, "opt",     ent_im_options},
  {BLTIN, "art",     ent_im_artifacts},
  {BLTIN, "iter",  ent_im_wand_iterate_images},
  {BLTIN, "read",     ent_im_read_image_to_wand},
  {BLTIN, "write",    ent_im_wand_write_image},
  {BLTIN, "disp",  ent_im_display_wand},
  {BLTIN, "exit",  ent_im_exit},
  {BLTIN, "method",  ent_im_wand_apply},
  {BLTIN, "join",  ent_im_wand_join},
  {BLTIN, "append",  ent_im_wand_append},
  {BLTIN, "combine",  ent_im_wand_combine},
  {BLTIN, "diff",  ent_im_wand_diff},
  {BLTIN, "deconstruct",  ent_im_wand_deconstruct},
  {BLTIN, "distort",  ent_im_wand_distort},
  {BLTIN, "f",  ent_im_wand_f},
  {BLTIN, "clut",  ent_im_wand_clut},
//   {BLTIN, "constitute",  ent_im_wand_image_constitute},
  {BLTIN, "pixels",  ent_im_wand_image_pixels},
  {BLTIN, "specs",  ent_im_wand_image_specs},
  {0, 0, 0},
};

