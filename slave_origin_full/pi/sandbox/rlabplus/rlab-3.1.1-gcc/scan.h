#ifndef SCAN_H
#define SCAN_H

/********************************************
scan.h
copyright 1991, Michael D. Brennan

This is a source file for mawk, an implementation of
the AWK programming language.

Mawk is distributed without warranty under the terms of
the GNU General Public License, version 2, 1991.
********************************************/

/*
 * Note: This source file has been modified extensively from
 * the original mawk version
 */

/*  Any parts not covered by the above are:
   Copyright (C) 1992, 1993  Ian R. Searle

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

#define  SC_SPACE               1
#define  SC_NL                  2
#define  SC_SEMI_COLON          3
#define  SC_FAKE_SEMI_COLON     4
#define  SC_LBRACE              5
#define  SC_RBRACE              6
#define  SC_QMARK               7
#define  SC_COLON               8
#define  SC_OR                  9
#define  SC_AND                10
#define  SC_PLUS               11
#define  SC_MINUS              12
#define  SC_MUL                13
#define  SC_DIV                14
#define  SC_MOD                15
#define  SC_POW                16
#define  SC_LPAREN             17
#define  SC_RPAREN             18
#define  SC_LBOX               19
#define  SC_RBOX               20
#define  SC_IDCHAR             21
#define  SC_DIGIT              22
#define  SC_DQUOTE             23
#define  SC_ESCAPE             24
#define  SC_COMMENT            25
#define  SC_EQUAL              26
#define  SC_NOT                27
#define  SC_LT                 28
#define  SC_GT                 29
#define  SC_COMMA              30
#define  SC_DOT                31
#define  SC_MATCH              32
#define  SC_DOLLAR             33
#define  SC_UNEXPECTED         34

#endif /* SCAN_H */
