/*\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/
\                                                                      /
/   Header file for internal      Copyright (c)  Dmitry A. Kazakov     \
\   pattern format                               St.Petersburg         /
/                                                Spring, 1993          \
\   (C, ANSI C)					                       /
/                                 Last revision :  17:50 30 Dec 1996   \
\                                                                      /
/   This  library  is  free software; you can redistribute it and/or   \
\   modify it under the terms of the GNU General Public  License  as   /
/   published by the Free Software Foundation; either version  2  of   \
\   the License, or (at your option) any later version. This library   /
/   is distributed in the hope that it will be useful,  but  WITHOUT   \
\   ANY   WARRANTY;   without   even   the   implied   warranty   of   /
/   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU   \
\   General  Public  License  for  more  details.  You  should  have   /
/   received  a  copy  of  the GNU General Public License along with   \
\   this library; if not, write to  the  Free  Software  Foundation,   /
/   Inc., 59 Temple Place - Suite 330, Boston, MA 02111-1307, USA.     \
\                                                                      /
/   As a special exception, if other files instantiate generics from   \
\   this unit, or you link this unit with other files to produce  an   /
/   executable, this unit does not by  itself  cause  the  resulting   \
\   executable to be covered by the GNU General Public License. This   /
/   exception  does not however invalidate any other reasons why the   \
\   executable file might be covered by the GNU Public License.        /
/                                                                      \
\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/*/

#ifndef	patran_h
#define	patran_h

/*>	Pattern internal representation

		Statement / Atom        Oct.	Mnemonic
		--------------------    ----	-----------------*/
#define		SUCCESS		        0000	/*	 S	 */
						/*>>Substrings<<<*/
#define		FIRST_ASCII_CHARACTER   0001	/*      soh	 */
#define		LAST_ASCII_CHARACTER    0177	/*      del	 */
						/*>>>Repeater<<<<*/
#define		MAX_REPEATER_COUNT       077	/*   small one	 */
#define		MAX_FINITE_REPEATER     0277	/*     0..63	 */
						/*>>>>Atoms<<<<<<*/
#define		END_OF_STRING	        0300	/*	 .	 */
#define		ANY		        0301	/*	 %	 */
#define		BLANK		        0302	/*	 +	 */
#define		DIGIT		        0303	/*	 #	 */
#define		UPPER_CASE_LETTER       0304	/*	 U	 */
#define		LOWER_CASE_LETTER       0305	/*	 W	 */
#define		LETTER		        0306	/*	 L	 */
#define		CHARACTER	        0307	/*	 C	 */
#define		FAILURE		        0310	/*	 F	 */
#define		END_OF_KEYWORD	        0311	/*	 _	 */
#define		NOOP		        0312	/* empty pattern */
#define		NEW_LINE	        0313	/*	 / 	 */
#define		CASE_DEAF_HOLERITH      0314	/*   <literal>	 */
#define		HOLERITH	        0315	/*   'literal'	 */
#define		ANY_OF			0316	/*   {literal}	 */
#define		MAX_HOLERITH	       	0400	/*  max length	 */
						/*>>Statements<<<*/
#define		FIRST_STATEMENT_CODE    0363	/*		 */
#define		ARB		        0363	/*	...	 */
#define		ASSIGN		        0364	/*       =	 */
#define		GO_TO		        0365	/*  goto label	 */
#define		FINITE_REPEATER	        0366	/*   0..2**32    */
#define		EXTERNAL_PATTERN        0367	/* call  pattern */
#define		QUERY		        0370	/*   no empty	 */
#define		DO_NOT_RETURN	        0371	/*	 :	 */
#define		INVERSE		        0372	/*      not	 */
#define		LITTLE_REPEATER	        0373	/*	 *	 */
#define		BIG_REPEATER	        0374	/*	 $	 */
#define		RIGHT_BRACE	        0375	/*	 )	 */
#define		OR		        0376	/*	 |	 */
#define		LEFT_BRACE	        0377	/*	 (       */
/*>

    In  most  cases  the statement or atom consists of exactly one byte.
    But there are several exceptions when a defined  sequence  of  bytes 
    follows  the  significant byte (one of described above). Here is the 
    list of the exceptions:

    (o)  The  CASE_DEAF_HOLERITH  or the HOLERITH bytes must be followed
         by  the  byte  counter (one byte) and the literal body with the 
	 length specified by that counter.

    (o)  The ANY_OF byte is followed by the pair of bytes. The first  of
         them  gives  5 high order bits of the lowest code of characters
         represented in the `any of' set. The second one gives  size  of
         the bit-map table that follows it. The table contains a bit set
         for each character from the `any of'.

    (o)  The  ASSIGN  byte  is followed by the four bytes containing the
         variable identifier number. The number format is machine depen-
         dent.

    (o)  The  GO_TO byte is followed by the four bytes containing signed
         relative offset to the target pattern. The offset is counted in
         the assumption that the next byte has offset 0. It has portable
         format. The low significant bit of the first byte is the offset
         sign.  Other 7 bits are low significant bits of the offset. The
         next byte contains the next 8 bits of the offset and so on.

    (o)  The EXTERNAL_PATTERN byte is followed by several bytes contain-
         ing the address of the target pattern. The number of bytes  and
         the address format are machine dependent.

    (o)  The FINITE_REPEATER byte is followed by the four bytes contain-
         ing the repetition count. The count format is the same  as  the
         offset format of GO_TO statement (see). It is portable.

<*/
#endif	/*>  patran_h  <*/
