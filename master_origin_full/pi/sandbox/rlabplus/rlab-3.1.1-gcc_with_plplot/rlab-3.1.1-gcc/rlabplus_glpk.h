// rlabplus extensions: glpk_simplex.c
extern Ent * ent_glpk_read_file (int nargs, Datum args[]);
extern Ent * ent_glpk_write_file (int nargs, Datum args[]);
extern Ent * ent_glpk_solve_lp  (int nargs, Datum args[]);
Bltin rlab_glpk_bltin[] = {
  // dloess
  {BLTIN, "read",  ent_glpk_read_file},
  {BLTIN, "write", ent_glpk_write_file},
  {BLTIN, "solve", ent_glpk_solve_lp},
  {0, 0, 0},
};
