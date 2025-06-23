//
// sys.op
//
// sys.op.zero
extern Ent *ent_Sys_Op_ZeroAbs (int nargs, Datum args[]);
extern Ent *ent_Sys_Op_ZeroRel (int nargs, Datum args[]);
Bltin rlab_sys_op_zero[] = {
{BLTIN, "abs", ent_Sys_Op_ZeroAbs},
{BLTIN, "rel", ent_Sys_Op_ZeroRel},
{0, 0, 0},
};

// sys.op
extern Ent *ent_Sys_Op_Max (int nargs, Datum args[]);
extern Ent *ent_Sys_Op_Min (int nargs, Datum args[]);
Bltin rlab_sys_op[] = {
  {BLTIN, "max", ent_Sys_Op_Max},
  {BLTIN, "min", ent_Sys_Op_Min},
  {0, 0, 0},
};


