/* code.h */

/*  This file is a part of RLaB ("Our"-LaB)
   Copyright (C) 1992, 1994  Ian R. Searle

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

#ifndef RLAB_CODE_H
#define RLAB_CODE_H

#include "rlab.h"
#include "list.h"
#include "symbol.h"

#include <stdio.h>

/*
 * The structure for a program
 */

struct _program
{
  int ncode;    /* The size of the progam array */
  Inst *prog;   /* The program array */
  Inst *progp;  /* Next free spot for code generation */
  int off;      /* Current instruction offset */
};

typedef struct _program Program;

extern void set_print_machine (int val);
extern int get_print_machine (void);
extern void set_line_nos (int val);
extern void set_use_pager (int val);
extern void set_code_pager (char *pager);
extern int get_progoff (void);
extern int set_progoff (void);

extern Program *program_Create (int n);
extern void program_Destroy (Program * p);
extern void program_Set (Program * p);
extern Inst * get_program_counter (void);

extern void init_machine (void);
extern void initcode (void);
extern void execute (Inst * p);
extern void execute_debug (Inst * p);
extern int find_lineno (void);
extern char *find_file_name (void);
extern Ent *get_eval_ret (void);

extern int code (int);
extern int codep (VPTR);
extern int coded (double);
extern void code_sp (int offset, int value);
extern int  get_code_sp (int offset);

extern void quit_code (void);
extern int delete_symbol_table (void);

extern void datum_stack_clean (void);
/* define op-codes */

#define STOP                  0
#define OP_PUSH_VAR           1
#define OP_PUSH_ARG           2
#define OP_PUSH_LOCAL_VAR     3
#define OP_PUSH_STATIC_VAR        82
#define OP_PULL_DATUM             83
#define OP_PUSH_DATUM             84
#define OP_ADD                4
#define OP_SUB                5
#define OP_MUL                6
#define OP_DIV                7
#define OP_LDIV               8
#define OP_NEGATE             9
#define OP_POWER             10
#define OP_ASSIGN            11
#define OP_FOR_LOOP_I        12
#define OP_FOR_LOOP_DONE     13
#define OP_FOR_LOOP_DONE_RET 14

#define OP_EL_MUL            15
#define OP_EL_DIV            16
#define OP_EL_LDIV           17
#define OP_EL_POWER          18
#define OP_PUSH_CONSTANT     19
#define OP_PUSH_iCONSTANT    20
#define OP_PRINT             21
#define OP_PRINT_ASSIGN           94

#define OP_GT                22
#define OP_LT                23
#define OP_EQ                24
#define OP_GE                25
#define OP_LE                26
#define OP_NE                27
#define OP_AND               28
#define OP_OR                29
#define OP_NOT               30
#define OP_ADDTO                  86
#define OP_SUBFROM                87
#define OP_EL_MUL_BY              88
#define OP_EL_DIV_BY              89
#define OP_IFSJMP            31
#define OP_FOR_LOOP          32
#define OP_FOR_THEN               85
#define OP_SWAP              33
#define OP_INC               34
#define OP_DEC               35
#define OP_POP               36
#define OP_VECTOR_CREATE     37
#define OP_VEC_APPEND        38
#define OP_MATRIX_VEC_SUB    39
#define OP_MATRIX_VEC_ASSIGN 40
#define OP_MATRIX_CREATE     41
#define OP_MATRIX_APPEND     42 
#define OP_MATRIX_ASSIGN     43
#define OP_MATRIX_SUB        44

#define OP_LIST_CREATE       45
#define OP_LIST_MEMB         46
#define OP_LIST_ASSIGN       47
#define OP_LIST_EL_CREATE    48
#define OP_FUNCTION_CALL     49
#define OP_FUNCTION_RETURN   50
#define OP_DEF_FUNC_RET      51
#define OP_BLTIN             52
#define OP_TRANSPOSE         53
#define OP_PUSH_STRING       54
#define OP_QUIT              55
#define OP_BREAK             56
#define OP_CONTINUE          57
#define OP_LINE_NO           58
#define OP_FILE_NAME         59

#define OP_FUNCTION_CALL_SELF 60
#define OP_JMP                61
#define OP_POP_CLEAN          62

#define OP_EMPTY_MATRIX_CREATE 63
#define OP_MATRIX_COL          64
#define OP_EL_TRANSPOSE        65

#define OP_RFILE               66
#define OP_RFILE_NAME          67

#define OP_HELP                68
#define OP_HELP_NAME           69

#define OP_PUSH_UNDEF          70
#define OP_EL_ADD              71
#define OP_EL_SUB              72

#define OP_SAVE_EVAL           73
#define OP_FUNCTION_CALL_1     74

#define OP_REQ_NAME            75

#define OP_OLIST_ASSIGN        76
#define OP_OLIST               77

#define OP_SYS_CMD             78

#define OP_PUSH_INTEGER        79

#define OP_EMPTY_MDE_CREATE    80
#define OP_FUNCTION_CALL_2          81
#define OP_ENDOFFILE                90
#define OP_JMP_ENDOFFILE            91

#define OP_DEF_CLASS_START          92
#define OP_DEF_CLASS_RET            93

#endif /* RLAB_CODE_H */
