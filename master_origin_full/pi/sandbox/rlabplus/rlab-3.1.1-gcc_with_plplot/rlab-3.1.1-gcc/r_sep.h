#ifndef RLABPLUS_SEP_H
#define RLABPLUS_SEP_H

#include "rlab.h"
#include "symbol.h"
#include "util.h"
#include "print.h"

extern Ent *ent_sep_background (int nargs, Datum args[]);
extern Ent *ent_sep_extract (int nargs, Datum args[]);

//
// sep
//
Bltin rlab_sep_bltin[] = {
  {BLTIN, "bkg", ent_sep_background},
  {BLTIN, "xtr", ent_sep_extract},
  {0, 0, 0},
};

#endif
