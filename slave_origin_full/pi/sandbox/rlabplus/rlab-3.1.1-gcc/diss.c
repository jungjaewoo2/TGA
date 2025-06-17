/* 
 * dissassem.c
 * Print human readable version of op-codes for inspection
 */

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

#include "rlab.h"
#include "code.h"
#include "bltin.h"
#include "print.h"
#include "util.h"

#include <stdio.h>

static Inst *pc;  /* Program counter (ptr) */
static int lineno;

/* **************************************************************
 * Functions to print various forms of op-codes.
 * ************************************************************** */
void print_op (Inst *p, char *s)
{
  fprintf (stderr, " %3d: %s\n", lineno++, s);
}

void print_code_ptr (Inst *p)
{
  /* pointers */
  fprintf (stderr, " %3d: %p\n", lineno++, (VPTR) (*p).ptr);
}

void print_code_string (Inst *p)
{
  /* pointers */
  fprintf (stderr, " %3d: %s\n", lineno++, (char *) ((*p).ptr));
}

void print_code_int (Inst *p)
{
  /* print ints */
  fprintf (stderr, " %3d: %d\n", lineno++, (*p).op_code);
}

void print_code_dval (Inst *p)
{
  /* doubles */
  fprintf (stderr, " %3d: %g\n", lineno++, (*p).d_val);
}

void print_code_var (Inst *p)
{
  /* print variables */
  if (((ListNode *) ((*p).ptr))->key != 0)
    fprintf (stderr, " %3d: %s\n", lineno++, ((ListNode *) ((*p).ptr))->key);
}

void diss_assemble (Inst * p, int pstop)
{
  int i;

  lineno = 1;
  pc = p;
  for (i = 0; i < pstop; i++)
  {
    switch ((*pc).op_code)
    {
    case OP_PUSH_VAR:
      print_op (pc++, "push var");
      print_code_var (pc++);
      i++;
      break;

    case OP_PUSH_ARG:
      print_op (pc++, "push arg var");
      print_code_int (pc++);
      i++;
      break;

    case OP_PUSH_LOCAL_VAR:
      print_op (pc++, "push local var");
      print_code_int (pc++);
      i++;
      break;

    case OP_ADD:
      print_op (pc++, "add");
      break;

    case OP_SUB:
      print_op (pc++, "sub");
      break;

    case OP_MUL:
      print_op (pc++, "multiply");
      break;

    case OP_DIV:
      print_op (pc++, "right divide");
      break;

    case OP_LDIV:
      print_op (pc++, "left divide");
      break;

    case OP_NEGATE:
      print_op (pc++, "negate");
      break;

    case OP_POWER:
      print_op (pc++, "power");
      break;

    case OP_ASSIGN:
      print_op (pc++, "assign");
      break;

    case OP_FOR_LOOP_I:
      print_op (pc++, "for-loop-init");
      print_code_int (pc++);
      i++;
      break;

    case OP_EL_MUL:
      print_op (pc++, "matrix-el-multiply");
      break;

    case OP_EL_DIV:
      print_op (pc++, "matrix-el-r-divide");
      break;

    case OP_EL_LDIV:
      print_op (pc++, "matrix-el-l-divide");
      break;

    case OP_EL_POWER:
      print_op (pc++, "matrix-el-power");
      break;

    case OP_PUSH_CONSTANT:
      print_op (pc++, "push constant");
      print_code_dval (pc++);
      i++;
      break;

    case OP_PUSH_iCONSTANT:
      print_op (pc++, "push imaginary constant");
      print_code_dval (pc++);
      i++;
      break;

    case OP_PRINT:
      print_op (pc++, "print");
      break;

    case OP_GT:
      print_op (pc++, "gt");
      break;

    case OP_LT:
      print_op (pc++, "lt");
      break;

    case OP_EQ:
      print_op (pc++, "eq");
      break;

    case OP_GE:
      print_op (pc++, "ge");
      break;

    case OP_LE:
      print_op (pc++, "le");
      break;

    case OP_NE:
      print_op (pc++, "ne");
      break;

    case OP_AND:
      print_op (pc++, "and");
      break;

    case OP_OR:
      print_op (pc++, "or");
      break;

    case OP_NOT:
      print_op (pc++, "not");
      break;

    case OP_IFSJMP:
      print_op (pc++, "ifsjmp");
      print_code_int (pc++);
      i++;
      break;

    case OP_FOR_LOOP:
      print_op (pc++, "for-loop");
      break;

    case OP_FOR_LOOP_DONE:
      print_op (pc++, "for-loop-done");
      break;

    case OP_FOR_LOOP_DONE_RET:
      print_op (pc++, "for-loop-done-ret");
      break;

    case OP_SWAP:
      print_op (pc++, "swap");
      break;

    case OP_INC:
      print_op (pc++, "inc");
      break;

    case OP_DEC:
      print_op (pc++, "dec");
      break;

    case OP_POP:
      print_op (pc++, "pop");
      break;

    case OP_POP_CLEAN:
      print_op (pc++, "pop clean");
      break;

    case OP_VECTOR_CREATE:
      print_op (pc++, "vector_create");
      print_code_int (pc++);
      i++;
      break;

    case OP_VEC_APPEND:
      print_op (pc++, "vector append");
      break;

    case OP_MATRIX_VEC_SUB:
      print_op (pc++, "matrix-vector sub");
      break;

    case OP_MATRIX_VEC_ASSIGN:
      print_op (pc++, "matrix-vector assign");
      break;

    case OP_MATRIX_CREATE:
      print_op (pc++, "matrix create");
      break;

    case OP_MATRIX_APPEND:
      print_op (pc++, "stack matrix");
      break;

    case OP_MATRIX_ASSIGN:
      print_op (pc++, "matrix assign");
      print_code_int (pc++);
      i++;
      break;

    case OP_MATRIX_SUB:
      print_op (pc++, "sub matrix");
      print_code_int (pc++);
      i++;
      break;

      case OP_LIST_CREATE:
        print_op (pc++, "create list");
        print_code_int (pc++);
        i++;
        break;

      case OP_LIST_MEMB:
        print_op (pc++, "list member");
        if (((*pc).op_code) == 1)
        {
          print_code_int (pc++);
          i++;
        }
        else
        {
          print_code_int (pc++);
        i++;
        print_code_string (pc++);
        i++;
      }
      break;

      case OP_LIST_ASSIGN:
        print_op (pc++, "list assign");
        if (((*pc).op_code) == 1)
        {
          print_code_int (pc++);
          i++;
        }
        else
      {
        print_code_int (pc++);
        i++;
        print_code_string (pc++);
        i++;
      }
      break;

    case OP_LIST_EL_CREATE:
      print_op (pc++, "list-el-create");
      print_code_ptr (pc++);
      i++;
      break;

    case OP_FUNCTION_CALL:
      print_op (pc++, "function call");
      print_code_int (pc++);
      i++;
      break;

    case OP_FUNCTION_CALL_1:
      print_op (pc++, "function call 1");
      print_code_int (pc++);
      i++;
      break;

    case OP_FUNCTION_CALL_SELF:
      print_op (pc++, "function call-self");
      print_code_int (pc++);
      i++;
      break;

    case OP_FUNCTION_RETURN:
      print_op (pc++, "function return");
      break;

    case OP_DEF_FUNC_RET:
      print_op (pc++, "default function return");
      break;

    case OP_TRANSPOSE:
      print_op (pc++, "matrix transpose");
      break;

    case OP_PUSH_STRING:
      print_op (pc++, "push_string");
      print_code_string (pc++);
      i++;
      break;

    case OP_QUIT:
      print_op (pc++, "quit");
      break;

    case OP_LINE_NO:
      fprintf (stderr, " %3d: %s", lineno++, "line # ");
      pc++;
      fprintf (stderr, " %d\n", (*pc++).op_code);
      i++;
      break;

    case OP_FILE_NAME:
      fprintf (stderr, " %3d: %s", lineno++, "file: ");
      pc++;
      fprintf (stderr, " %s\n", (char *) (*pc++).ptr);
      i++;
      break;

    case OP_JMP:
      print_op (pc++, "jmp");
      print_code_int (pc++);
      i++;
      break;

    case OP_EMPTY_MATRIX_CREATE:
      print_op (pc++, "create_empty_matrix");
      break;

    case OP_MATRIX_COL:
      print_op (pc++, "matrix_reshape_col");
      break;

    case OP_EL_TRANSPOSE:
      print_op (pc++, "matrix el-transpose");
      break;

    case OP_RFILE:
      print_op (pc++, "rfile command");
      break;

    case OP_RFILE_NAME:
      print_op (pc++, "rfile-name command");
      print_code_string (pc++);
      i++;
      break;

    case OP_REQ_NAME:
      print_op (pc++, "require command");
      print_code_string (pc++);
      i++;
      break;

    case OP_HELP:
      print_op (pc++, "help command");
      break;

    case OP_HELP_NAME:
      print_op (pc++, "help-name command");
      print_code_string (pc++);
      i++;
      break;

    case OP_PUSH_UNDEF:
      print_op (pc++, "push-UNDEF");
      break;

    case OP_EL_ADD:
      print_op (pc++, "el-add");
      break;

    case OP_EL_SUB:
      print_op (pc++, "el-sub");
      break;

    case STOP:
      print_op (pc++, "stop");
      break;

    case OP_SAVE_EVAL:
      print_op (pc++, "save_eval");
      break;

    case OP_OLIST_ASSIGN:
      print_op (pc++, "olist_assign");
      break;

    case OP_OLIST:
      print_op (pc++, "olist");
      print_code_int (pc++);
      i++;
      break;

    case OP_SYS_CMD:
      print_op (pc++, "sys_cmd");
      print_code_string (pc++);
      i++;
      break;

    default:
      fprintf (stderr, "Invalid op-code: %d\n", (*pc).op_code);
      fflush (stderr);
      rerror ("invalid op-code for diss-assembler");
      break;
    }
  }
}
