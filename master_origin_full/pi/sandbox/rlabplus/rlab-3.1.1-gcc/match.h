/*
\                                                                      /
/   Header file for pattern         Copyright (c)  Dmitry A. Kazakov   \
\   matching (externals)                           St.Petersburg       /
/                                                  Spring, 1993        \
\                                                                      /
/   (C, ANSI C, C++)              Last revision :  15:59 22 Jan 2000   \
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

  \                          /
   /	  USER INTERFACE      \
   \	External variables    /

   ----------
   UserAlphas	- list of user defined alphas

   This  is a NUL terminated string that contains the list of symbols to
   be interpreted as valid characters  for  identifiers  of  labels  and
   variables. Normally UserAlphas contains 0. It means that only letters
   and digits are valid. If for example UserAlphas = "_#" then both  un-
   derline  and number characters are valid too. UserAlphas's characters
   are assumed to be letters. Notice the difference between letters  and
   digits: an indentifier name always starts with a letter.

   --------
   UserBase	- base of user digits

   This is an integer that contains the current base for the DIGIT atom.
   Normally it is set to 10, so DIGIT means one of the following charac-
   ters 0123456789. Setting it to, say, 2 causes DIGIT be 01. The  User-
   Base  could  be greater than 10. For example, setting it to 16 causes
   DIGIT be 0123456789abcdefABCDEF. Note, that  the  CHARACTER  and  the
   BREAK atoms do not depend on UserBase value.

   ----------
   MatchError	- last error code

   This variable is set after each call of match,  patran  or  patmaker.
   It contains the same value as one returned by match and patran.


   USER INTERFACE
   External functions

   All mentioned functions are specified by pointers. It is possible  to
   reset such pointer to the appropriate code or to set it to NULL.

   -----------
   GetNextLine	- Called to pick up next line to match

      Line	- address of the next line address
      Length	- length of the next line (must be >= 0)

   Returns :

      0 if there is no a line (end  of  file  equivalent).  Other  value
      means that Line and Length were correctly set.

   ---------------
   GetPreviousLine  - Called to return back to the previous line

      Line	- address of the previous line address
      Length	- length of the previous line

   Returns :

      0 if the previous line cannot be returned. In this  case  matching
      process  will  be terminated with an error. Other value means that
      Line and Length were correctly set.

   ------------------
   GetExternalPattern  - Find for an external pattern

      Name	- address of a pattern name
      Length	- name length

   Returns :

      Pattern or NULL.

   This  function is called when the pattern translator discovers an un-
   defined identifier. This identifier is passed  to  the  GetExt...  to
   provide an embedded pattern feature. If the GetExt... returns NULL it
   means that there is no such pattern and the pattern translation  will
   fail.

   -------------
   GetVariableId  - Find for a variable

      Name	- address of a variable name
      Length	- name length

   Returns :

      VariableId or 0.

   This  function  is  called  when  the pattern translator discovers an
   equation. The user defined function must return the variable id - any
   non-zero  number  that  will be passed to the AssignVariable (see). 0
   means that there is no such variable and the translation will fail.

   --------------
   AssignVariable  - Assign value on successful matching

      VariableId	- id (returned by GetVariableId)
      Offset		- offset to beginning of matched part of string
      Length		- length of matched part

   This  function  is  called each time when a pattern after an equation
   sign is matched. The VariableId is the same as one  returned  by  the
   GetVariableId.  the  Offset  and  the Length define which part of the
   string was matched. Note that the string  may  occupy  several  lines
   (see  the  GetNextLine) so the Offset means an offset from the begin-
   ning of the first line including line ends (a line end is encountered
   as one character).

   ----------------
   DeAssignVariable  - Deassign value on failure

      VariableId	- id (returned by GetVariableId)

   A  call  to the AssignVariable does not mean that a pattern after the
   equation sign was finally matched. Future matching might fail and the
   previously assigned value must be discarded. For this purpose the De-
   AssignVariable is called. The VariableId is the same as one  returned
   by the GetVariableId.

   ---------------
   AllocatePattern  - Memory allocator for the pattern body

      Size		- required size in bytes

   Returns :

      Address of the memory block (char *) or 0

   The  function is called by the patinit (see below) to allocate memory
   for the pattern being translated. AllocatePattern must return  either
   address of the allocated memory or 0 indicating failure.

<*/
#ifndef	match_h
#define	match_h

#ifdef __cplusplus

extern "C"		// This allows you to use a C++ compiler
{

#endif  /*> C++ <*/

#ifndef NON_ANSI
extern int   	     (* GetNextLine)        (const char ** Line, int * Length);
extern int   	     (* GetPreviousLine)    (const char ** Line, int * Length);
extern char *	     (* AllocatePattern)    (unsigned long Size);
extern char *	     (* GetExternalPattern) (char *        Name, int   Length);
// extern unsigned long (* GetVariableId)	    (char *        Name, int   Length);
/*extern void          (* AssignVariable)     (
					       unsigned long VariableId,
					       unsigned long Offset,
					       unsigned long Length
					    );*/
// extern void          (* DeAssignVariable)   (unsigned long VariableId);
#else
extern int   	     (* GetNextLine)        ();
extern int   	     (* GetPreviousLine)    ();
extern char *	     (* AllocatePattern)    ();
extern char *	     (* GetExternalPattern) ();
extern unsigned long (* GetVariableId)	    ();
extern void          (* AssignVariable)     ();
extern void          (* DeAssignVariable)   ();
#endif	/*> NON_ANSI <*/

extern char *	     UserAlphas;
extern int	     UserBase;
extern int	     MatchError;

/*>	Status codes	<*/

#define 	MATCH_SUCCESS			  1
#define		MATCH_FAILURE			  0
#define		MATCH_WRONG_PATTERN_FORMAT	 -1
#define		MATCH_NO_DYNAMIC_MEMORY		 -2
#define 	MATCH_BRACE_ERROR		 -3
#define		MATCH_MISSING_QUOTATION		 -4
#define		MATCH_TOO_BIG_REPEATER		 -5
#define		MATCH_DUPLICATE_LABEL		 -6
#define 	MATCH_UNRECOGNIZED_CHARACTER	 -7
#define		MATCH_MISSING_RIGHT_BRACE	 -8
#define		MATCH_UNRECOGNIZED_KEYWORD	 -9
#define		MATCH_TOO_LARGE_PATTERN		-10
#define	 	MATCH_TOO_LARGE_STRING		-11
#define	 	MATCH_CANNOT_RETURN		-12
#define		MATCH_UNDEFINED_VARIABLE	-13
#define		MATCH_POSSIBLE_INDEFINITE_LOOP	-14
#define		MATCH_RESERVED_KEYWORD		-15
#define		MATCH_INTERNAL_ERROR		-16

/*>	External functions	<*/

/*>

   match -- match a pattern

	Length	- of a character string to be matched
	String	- points to the string (maybe not NUL terminated)
	Pointer - number of 1st character to be matched (0..Length)
	Pattern	- points to a pattern (in internal format)

   Returns :
		     MATCH_SUCCESS - Success
		     MATCH_FAILURE - Failure of matching
	MATCH_WRONG_PATTERN_FORMAT - Wrong pattern format
	     MATCH_TOO_LONG_STRING - String length exceeds 1Gbyte
	       MATCH_CANNOT_RETURN - Cannot return to previously left string
	   MATCH_NO_DYNAMIC_MEMORY - Dynamic memory overflow
	      MATCH_INTERNAL_ERROR - Unrecognized internal error

   After successful completion String [*Pointer] will be the 1st charac-
   ter  following matched substring. Note that after successful matching
   the new line atom (/) Pointer will be still counted as an offset from
   beginning  of  the String. The offset includes line ends: +1 for each
   end.

<*/
#ifndef NON_ANSI
extern int match
(
   int	      	Length,
   const char *	String,
   int        * Pointer,
   const char *	Pattern
);
#else
extern int match ();
#endif	/*> NON_ANSI <*/
/*>

   patinit -- pattern constructor (see also patmaker)

	Length	- of a character string containing a pattern to be translated
	String	- points to the string (maybe not NUL terminated)
	Pointer - number of starting character (0..Length)

   Returns :

	The address of the translated pattern or 0

   After  successful  completion  String  [*Pointer]  will  be  the  1st
   character following translated pattern. The AllocatePattern  function
   is  called  to  allocate  translated  pattern. On success the pattern
   address  is  returned.  Otherwise  the  result  is  0. The MatchError
   variable can be consulted for the actual error code.

<*/
#ifndef NON_ANSI
extern const char * patinit
(
   int		Length,
   const char *	String,
   int        *	Pointer
);
#else
extern char * patinit ();
#endif	/*> NON_ANSI <*/
/*>

   patran -- pattern translator

	Length	- of a character string containing a pattern to be translated
	String	- points to the string (maybe not NUL terminated)
	Pointer - number of starting character (0..Length)
	Pattern	- address of a buffer for translated pattern
	Size	- its size

   Returns :
		     MATCH_SUCCESS - Success
		 MATCH_BRACE_ERROR - Right brace does not match the left one
	     MATCH_DUPLICATE_LABEL - Duplicate definition of a label
		     MATCH_FAILURE - Failure. There is no any pattern
	   MATCH_MISSING_QUOTATION - Missing closing quotation marks
	 MATCH_MISSING_RIGHT_BRACE - One or more missing right braces
	   MATCH_NO_DYNAMIC_MEMORY - Dynamic memory overflow
    MATCH_POSSIBLE_INDEFINITE_LOOP - Repetition of a pattern that matches null
	    MATCH_TOO_BIG_REPEATER - Repetition count is too big
	   MATCH_TOO_LARGE_PATTERN - Pattern is too large for provided buffer
	    MATCH_TOO_LARGE_STRING - Pattern length exceeds 1Gbyte
	  MATCH_UNDEFINED_VARIABLE - Immediate assignment to unknown variable
      MATCH_UNRECOGNIZED_CHARACTER - Unrecognized character
	MATCH_UNRECOGNIZED_KEYWORD - Unrecognized keyword or undefined label
	    MATCH_RESERVED_KEYWORD - Empty or conflict variable name

   After successful completion String [*Pointer] will be the 1st charac-
   ter  following  translated  pattern. Size will contain actual size of
   the pattern i.e. number of used out bytes of the Pattern buffer.

<*/
#ifndef NON_ANSI
extern int patran
(
   int		Length,
   const char *	String,
   int        *	Pointer,
   char       *	Pattern,
   int        *	Size
);
#else
extern int patran ();
#endif	/*> NON_ANSI <*/
/*>

   patmaker -- pattern constructor

	String	- points to the NUL terminated string

   Returns :

	Pointer to translated pattern or NULL

   The String is translated and a memory block is requested to  allocate
   the  pattern  body.  Any syntax or other error leads to the NULL as a
   result. Note that whole String must be recognized as a pattern. Error
   code  is returned to the MatchError variable.  You should member that
   error may be resulted from malloc's failure.

<*/
#ifndef NON_ANSI
extern const char * patmaker (const char * String);
#else
extern char * patmaker ();
#endif	/*> NON_ANSI <*/
/*>

   freedm -- releases dynamic memory occupied by ASS and MSS

   After execution the freedm ASS and MSS pointers will be  restored  to
   previous  values  (First_ASS_ptr,  First_ASS_used,  ...).  All  stack
   frames allocated later will be released.

			      Will be released >>>>>>>>>>>>>>>>>>>>>>>
		._________.      ._________.      ._________.
		|         |      |         |      |         |
	    --->|  frame  |----->|  frame  |----->|  frame  |----> nil
		|_________|      |_________|      |_________|
		       |              |
       First_xSS_ptr __|    xSS_ptr __|                     x = A or M

<*/
#ifndef NON_ANSI
extern void freedm (void);
#else
extern void freedm ();
#endif	/*> NON_ANSI <*/

#ifdef	__cplusplus
}
#endif	/*> C++ <*/

#endif	/*>  match_h  <*/
