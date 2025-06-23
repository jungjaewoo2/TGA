// sparskit.h

//
// Defines and declarations for a C => fortran, RLaB => Sparskit interface.
//
// by Marijan Kostrun, (c) VI 2005-2013

#ifndef RLAB_SPARSKIT_H
#define RLAB_SPARSKIT_H

#include "fi.h"
#include "sparse.h"

#ifdef HAVE_FORTRAN_UND_BACK
// INOUT Module
#define OPENFILE openfile_
#define CLOSFILE closfile_
#define PRTMT prtmt_
#define SMMS smms_
#define DUMP dump_
#define PSPLTM pspltm_
// ILUT Module
#define PGMRES pgmres_
#define ILUT ilut_
#define ILUTP ilutp_
#define ILUD ilud_
#define ILUDP iludp_
#define ILUK iluk_
// ITERS Module
#define CG cg_
#define CGNR cgnr_
#define BCG bcg_
#define BCGSTAB bcgstab_
#define TFQMR tfqmr_
#define GMRES gmres_
#define FGMRES fgmres_
#define DQGMRES dqgmres_
#define DBCG dbcg_
#define DISTDOT distdot_
#define RUNRC runrc_
// MATVEC Module
#define AMUX amux_
#define ATMUX atmux_
#define SPLUSOL splusol_
#define LUTSOL lutsol_
#endif

#ifdef HAVE_FORTRAN_UND_FRONT
// INOUT Module
#define OPENFILE _openfile
#define CLOSFILE _closfile
#define PRTMT _prtmt
#define SMMS _smms
#define DUMP _dump
#define PSPLTM _pspltm
// ILUT Module
#define PGMRES _pgmres
#define ILUT _ilut
#define ILUD _ilud
#define ILUTP _ilutp
#define ILUDP _iludp
#define ILUK _iluk
// ITERS Module
#define CG _cg
#define CGNR _cgnr
#define BCG _bcg
#define BCGSTAB _bcgstab
#define TFQMR _tfqmr
#define GMRES _gmres
#define FGMRES _fgmres
#define DQGMRES _dqgmres
#define DBCG _dbcg
#define DISTDOT _distdot
#define RUNRC _runrc
// MATVEC Module
#define AMUX _amux
#define ATMUX _atmux
#define SPLUSOL _splusol
#define LUTSOL _lutsol
#endif

#ifdef HAVE_FORTRAN_UPPERCASE
// Do nothing, the existing code is OK
#endif

#ifdef HAVE_FORTRAN_LOWERCASE
// INOUT Module
#define OPENFILE openfile
#define CLOSFILE closfile
#define PRTMT prtmt
#define SMMS smms
#define DUMP dump
#define PSPLTM pspltm
// ILUT Module
#define PGMRES pgmres
#define ILUT ilut
#define ILUD ilud
#define ILUTP ilutp
#define ILUDP iludp
#define ILUK iluk
// ITERS Module
#define CG cg
#define CGNR cgnr
#define BCG bcg
#define BCGSTAB bcgstab
#define TFQMR tfqmr
#define GMRES gmres
#define FGMRES fgmres
#define DQGMRES dqgmres
#define DBCG dbcg
#define DISTDOT distdot
#define RUNRC runrc
// MATVEC Module
#define AMUX amux
#define ATMUX atmux
#define SPLUSOL splusol
#define LUTSOL lutsol
#endif


//
// INOUT Module
//
extern int OPENFILE (char *, int *, int *);
extern int CLOSFILE (int *);
extern int PRTMT (int *, int *, double *, int *, int *, char *, char *, char *, char *, char *,
                  int *, int *, int *);
extern int SMMS (int *, int *, int *, int *, double *, int *, int *, int *);
extern int DUMP (int *, int *, int *, double *, int *, int *, int *);
extern int PSPLTM (int *, int *, int *, int *, int *, char *, int *, double *,
                   char *, int *, int *, int *);
//
// ILUT Module
//
extern int PGMRES();
extern int ILUT(int *, double *, int *, int *, int *, double *, double *, int *, int *, int *, double *,
                int *, int *);    // 1
extern int ILUTP(int *, double *, int *, int *, int *, double *, double *, int *, int *, int *,
                 double *, int *, int *);   // 2
extern int ILUD(int *, double *, int *, int *, int *, double *, double *, int *, int *, int *, double *,
                int *, int *);    // 3
extern int ILUDP(int *, double *, int *, int *, int *, double *, double *, int *, int *, int *,
                 double *, int *, int *);   // 4
extern int ILUK(int *, double *, int *, int *, int *, double *, double *, int *, int *, int *, double *,
                int *, int *);    // 5
//
// ITERS Module
//
extern int CG();      // 1
extern int CGNR();    // 2
extern int BCG();     // 3
extern int BCGSTAB(); // 4
extern int TFQMR();   // 5
extern int GMRES();   // 6
extern int FGMRES();  // 7
extern int DQGMRES(); // 8
extern int DBCG();    // 9
extern int DISTDOT();
extern int RUNRC(int *, double *, double *, int *, double *, double *, double *, int *, int *, double *,
                 int *, int *, int *);
//
// MATVEC module
//
extern int AMUX();
extern int ATMUX();
extern int SPLUSOL();
extern int LUTSOL();


#endif  /* RLAB_SPARSKIT_H */
