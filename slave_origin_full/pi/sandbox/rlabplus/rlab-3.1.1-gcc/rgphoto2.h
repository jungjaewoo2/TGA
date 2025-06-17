/* rgphoto2.h */

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

extern Ent *ent_gp_camera_init    (int nargs, Datum args[]);
extern Ent *ent_gp_camera_exit    (int nargs, Datum args[]);
extern Ent *ent_gp_camera_config  (int nargs, Datum args[]);
extern Ent *ent_gp_camera_info    (int nargs, Datum args[]);
extern Ent *ent_gp_camera_capture (int nargs, Datum args[]);
extern Ent *ent_gp_camera_config_options (int nargs, Datum args[]);
extern Ent *ent_gp_info_supported_cameras (int nargs, Datum args[]);

Bltin rlab_gphoto2_bltin[] = {
  {BLTIN, "init",       ent_gp_camera_init},
  {BLTIN, "exit",       ent_gp_camera_exit},
  {BLTIN, "config",     ent_gp_camera_config},
  {BLTIN, "config_options", ent_gp_camera_config_options},
  {BLTIN, "info",       ent_gp_camera_info},
  {BLTIN, "capture",       ent_gp_camera_capture},
  {BLTIN, "supported",  ent_gp_info_supported_cameras},
  {0, 0, 0},
};

