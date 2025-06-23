/*\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/*\
\                                                                      /
/   Header file for pattern         Copyright (c)  Dmitry A. Kazakov   \
\   type                                           St.Petersburg       /
/                                                  Spring 1993         \
\   (C++)                                                              /
/                                 Last revision :  14:13 20 Oct 2001   \
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
\*\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/*/

#ifndef	pattern_h
#define	pattern_h
//
// This header file is provided for C++. It defines the type pattern and
// some  useful  operators  on it. All of them are implemented as inline
// functions because all the work is made in the match.c. 
//
//			M  E  T  H  O  D  S
//
// Kind of      Result        Notation	   Arguments
// -----------  ------------  -----------  -----------------------------
// friend       int	      operator |=  const char *, const pattern&
// friend       int	      operator >>  const char *, const pattern&
// constructor	       	      pattern	   const char *
// destructor	      	     ~pattern 	
// member       const char *  Pattern ()		
// -----------  ------------  -----------  -----------------------------
//
//			C  O  M  M  E  N  T  S
//
// The constructor  allocates  the  buffer  for  the  pattern  body  and
// translates  it  from  the  textual  form  to the internal format. Any
// translation error lead to setting the field Body to 0.  That  pattern
// will match nothing. The variable MatchError  should  be  observed  to
// determine the actual error code (see match.h for more information). 
//
// The pattern body can  be  directly  accessed  as  Pattern  (X).  This
// function returns the pointer to the first character of the  pattern's
// internal representation. 
//
// The  operator  <string>  |=  <pattern> matches the specified <string>
// against the <pattern> and returns 1  if  the  <pattern>  matches  the
// whole <string>.
//
// The operator <string> >> <pattern> matches a prefix of the  specified
// <string> and returns the number of matched characters.  So,  you  can
// continue  matching  process  with  the  aid  of  another pattern. For
// example : 
//
//	int  Index;
//             . . .
//	Index  = 0;                          // From the beginning
//	Index += &Line [Index] >> Patern1;   // Skip by Pattern1
//	Index += &Line [Index] >> Patern2;   // Skip by Pattern2
//	     . . .
//
// The  destructor  releases  the  heap  memory allocated by the pattern
// body.
// 
// It is possible but not desirable to match  strings  with  dynamically
// created  patterns,  like  String  >> "letter $character:". Because it
// leads  to  pattern  translation  each  time  when  the  statement  is
// executed. String to pattern conversion and pattern constructor  allow
// declarations of pattern objects:
//
//	const pattern	Identifier ("letter $character:");
//
#include	<stdlib.h>		// For malloc and free
#include	<string.h>  		// For strlen
#include	"match.h"  		// For match and patmaker
//
// The  definition  of  the  class pattern. The only member is protected
// field  Body.  It  contains  the  pointer  to  the  internal   pattern
// representation.  The  pattern body is allocated when an object of the
// pattern class is created and released when the object is removed. 
//
class pattern
{
public :
//
// Constructor -- From string
//
//	String	- To be translated
//
   pattern (const char * String) : Body (patmaker (String)) {}
//
// operator |= -- Match the whole string against the pattern
//
//	String	- Pointer to the string to be matched
//	Pattern	- Pattern
//
// Returns :
//
//	[0]  The pattern does not match the whole string
//	[1]  The pattern matches the whole string
//
   friend int operator |=
   (
      const char     *	String,
      const pattern&	Pattern
   );
//
// operator >> -- Match a part of the string against the pattern
//
//	String	- Pointer to the string to be matched
//	Pattern	- Pattern
//
// Returns :
//
//	Number of matched string bytes
// 
   friend int operator >>
   (
      const char     *	String,
      const pattern&	Pattern
   );
//
// Destructor
//
  ~pattern ()
   {
      if (Body) free ((void *) Body);
   }
//
// Pattern -- The internal pattern representation
//
// Returns :
//
//	Pointer to the string representing the pattern
//
   const char * Pattern () const
   {
      return Body;
   }

private :
   const char *	Body;  	// Pointer to the internal representation
//
// Preventing automatic copying of a pattern ...
//
   pattern (const pattern&);		// Declared but not defined 
   void operator = (const pattern&);  	// Declared but not defined

}; // pattern

inline int operator |= (const char * String, const pattern& Pattern)
{
   int	Index  = 0;
   int	Length = strlen (String);

   return
   (  0 != Pattern.Pattern ()
   && MATCH_SUCCESS == match (Length, String, &Index, Pattern.Body)
   && Length == Index
   );
};

inline int operator >> (const char * String, const pattern& Pattern)
{  // Matching of string prefix returns number of matched characters
   int	Index = 0;

   if (0 != Pattern.Pattern ())
   {  // Here is a body
      (void) match (strlen (String), String, &Index, Pattern.Body);
   }
   return (Index);
};

#endif	/*>  pattern_h  <*/
