//
// python for rlabplus
//

extern Ent *ent_Py_Initialize (int nargs, Datum args[]);
extern Ent *ent_Py_GetVersion (int nargs, Datum args[]);
extern Ent *ent_Py_Finalize (int nargs, Datum args[]);
extern Ent *ent_PyRun_SimpleString (int nargs, Datum args[]);
extern Ent *ent_PyRun_File (int nargs, Datum args[]);
extern Ent *ent_PyRun_EvaluateExpression (int nargs, Datum args[]);
extern Ent *ent_Py_ManipulateObject (int nargs, Datum args[]);

Bltin rlab_python_bltin[] = {
  {BLTIN, "init", ent_Py_Initialize},
  {BLTIN, "ver",  ent_Py_GetVersion},
  {BLTIN, "kill", ent_Py_Finalize},
  {BLTIN, "cmd",  ent_PyRun_SimpleString},
  {BLTIN, "file", ent_PyRun_File},
  {BLTIN, "eval", ent_PyRun_EvaluateExpression},
  {BLTIN, "var", ent_Py_ManipulateObject},
  {0, 0, 0},
};
