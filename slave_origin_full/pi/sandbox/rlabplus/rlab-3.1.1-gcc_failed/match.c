//  Pattern matching                Copyright (c)  Dmitry A. Kazakov
//  freedm, match, patinit,                        St.Petersburg
//  patran, patmaker                               Autumn, 1993
//
//  (C, ANSI C)                   Last revision :  17:49 30 Dec 1996
//
//  This  library  is  free software; you can redistribute it and/or
//  modify it under the terms of the GNU General Public  License  as
//  published by the Free Software Foundation; either version  2  of
//  the License, or (at your option) any later version. This library
//  is distributed in the hope that it will be useful,  but  WITHOUT
//  ANY   WARRANTY;   without   even   the   implied   warranty   of
//  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU
//  General  Public  License  for  more  details.  You  should  have
//  received  a  copy  of  the GNU General Public License along with
//  this library; if not, write to  the  Free  Software  Foundation,
//  Inc., 59 Temple Place - Suite 330, Boston, MA 02111-1307, USA.
//
//  As a special exception, if other files instantiate generics from
//  this unit, or you link this unit with other files to produce  an
//  executable, this unit does not by  itself  cause  the  resulting
//  executable to be covered by the GNU General Public License. This
//  exception  does not however invalidate any other reasons why the
//  executable file might be covered by the GNU Public License.       

//
//
//   Symbol  DEBUG  defines  debugging level i.e. how much verbose must be
//   the program. A working version should be compiled without  DEBUG  de-
//   fined.
//
//        #define DEBUG 9
//
//   Symbol  PTR_SIZE forces exactly specified number bytes to be used for
//   internal data pointers representation.  By default the number of used
//   bytes is either 4 (for 16/32-bit machines) or 8 (for more than 32-bit
//   machines).  You can define it as 2  for  a 16-bit  processor  without
//   memory control unit.
//
//        #define PTR_SIZE 2
//
// #define DEBUG 9

#ifdef DEBUG
#include        <assert.h>
#endif
#ifndef NON_ANSI
// #include        <gc.h>      /* For malloc                           */
#include        <string.h>      /* For strlen, strchr                   */
#else
#define const                   /* Ignore const keyword                 */
#endif /*> NON_ANSI <*/

#include        <ctype.h>       /* For isalpha and so on                */
#include        <string.h>      /* For strlen, strchr                   */
#include        <stdio.h>       /* Just to define NULL                  */

#include	"sizeof.h"	/* Machine properties			*/

#include        "match.h"       /* Definitions for myself               */
#include        "patran.h"      /* Internal representation of patterns  */
/*>

   Definition of useful types

<*/
#ifndef PTR_SIZE
#if MachineWord == 8
#define         PTR_SIZE        8
#else
#define         PTR_SIZE        4
#endif	/*> MachineWord == 8 <*/
#endif  /*> PTR_SIZE <*/

typedef         char                    boolean;
typedef         unsigned char           byte;
typedef         unsigned short          halfword;
typedef         unsigned int            natural;
#if MachineWord == 2
typedef         unsigned long           word;
#else
                                /* 32 or more bit machine		*/
typedef         unsigned int            word;
#endif	/*> MachineWord == 2 <*/

#define         TRUE                    1
#define         FALSE                   0
#define         WORD_SIZE               4

#define         APOSTROPHE              '\''
#define         CIRCUMFLEX              '^'
#define         EQUATION                '='
#define         LEFT_ANGLE              '<'
#define         NUL                     '\000'
#define         RABBIT_EARS             '"'
#define         RIGHT_ANGLE             '>'
#define         SPACE                   ' '
#define         TAB                     '\t'
#define         LF                      '\n'

/*>     Stack control parameters        <*/

#if MachineWord == 2
#define         MAX_STACK_FRAME         512
#else
#define         MAX_STACK_FRAME         4096
#endif

#define         STACK_FRAME_SIZE        MAX_STACK_FRAME + PTR_SIZE * 2

/*>     Some masks for managing statement header tag    <*/

#define         HEADER_TEST_MASK        0x81
#define         HEADER_SLICE_MASK      ~0x81

/*>     Other constants         <*/

#define         MAX_INT_DIV_10          214748364L
#define         LAST_MSS_MASK           128
#define         KEY_LIST_SIZE           38
#define         MAX_STRING_LENGTH       0x4FFFFFFFL

/*^v^v^v^v^v^v^v^v^v^v^v^v^v^v^ MACROS ^v^v^v^v^v^v^v^v^v^v^v^v^v^v^v^v^*/
/*>

   PopPtr  -- Pop pointer in reverse order
   PushPtr -- Push pointer in reverse order

        BPtr    - Pointer to byte following the last byte

   This macro copies subsequent bytes in reverse order. After completion
   pointer  will  point to  the first byte.  The copy is  stored  in the
   Take.Byte array.

<*/
#if PTR_SIZE == 2
#define         PopPtr(BPtr)						\
                {                                                       \
                   Take.Byte [1] = *--BPtr;                             \
                   Take.Byte [0] = *--BPtr;                             \
                }
#define         PushPtr(BPtr, AVal)					\
                {                                                       \
		   Take.Ptr = AVal;					\
                   *--BPtr = Take.Byte [1];                             \
                   *--BPtr = Take.Byte [0];                             \
                }
#else
#if PTR_SIZE == 4
#define         PopPtr(BPtr)						\
                {                                                       \
                   Take.Byte [3] = *--BPtr;                             \
                   Take.Byte [2] = *--BPtr;                             \
                   Take.Byte [1] = *--BPtr;                             \
                   Take.Byte [0] = *--BPtr;                             \
                }
#define         PushPtr(BPtr, AVal)					\
                {                                                       \
		   Take.Ptr = AVal;					\
                   *--BPtr = Take.Byte [3];                             \
                   *--BPtr = Take.Byte [2];                             \
                   *--BPtr = Take.Byte [1];                             \
                   *--BPtr = Take.Byte [0];                             \
                }
#else
#define         PopPtr(BPtr)						\
                {                                                       \
                   Take.Byte [7] = *--BPtr;                             \
                   Take.Byte [6] = *--BPtr;                             \
                   Take.Byte [5] = *--BPtr;                             \
                   Take.Byte [4] = *--BPtr;                             \
                   Take.Byte [3] = *--BPtr;                             \
                   Take.Byte [2] = *--BPtr;                             \
                   Take.Byte [1] = *--BPtr;                             \
                   Take.Byte [0] = *--BPtr;                             \
                }
#define         PushPtr(BPtr, AVal)					\
                {                                                       \
		   Take.Ptr = AVal;					\
                   *--BPtr = Take.Byte [7];                             \
                   *--BPtr = Take.Byte [6];                             \
                   *--BPtr = Take.Byte [5];                             \
                   *--BPtr = Take.Byte [4];                             \
                   *--BPtr = Take.Byte [3];                             \
                   *--BPtr = Take.Byte [2];                             \
                   *--BPtr = Take.Byte [1];                             \
                   *--BPtr = Take.Byte [0];                             \
                }
#endif
#endif
/*>

   PickUpWord -- Assembly a word from four subsequential bytes

        BPtr    - Pointer to the fisrt byte

   This  macro  takes  four  subsequential  bytes  and  stores them into
   Take.Byte variable. Pointer is advanced to  the  byte  following  the
   last byte.

<*/
#define         PickUpWord(BPtr)                                        \
                (                                                       \
                   Take.Byte [0] = *BPtr++,                             \
                   Take.Byte [1] = *BPtr++,                             \
                   Take.Byte [2] = *BPtr++,                             \
                   Take.Byte [3] = *BPtr++                              \
                )
/*>

   PushWord -- Put four bytes

        BPtr    - Pointer to the fisrt destination byte

   This  macro  stores Take.Byte into four subsequent bytes. The pointer
   is advanced to the byte following the last destination byte.

<*/
#define         PushWord(BPtr)                                          \
                {                                                       \
                   *BPtr++ = Take.Byte [0];                             \
                   *BPtr++ = Take.Byte [1];                             \
                   *BPtr++ = Take.Byte [2];                             \
                   *BPtr++ = Take.Byte [3];                             \
                }
/*>

   PutWord -- Store a word into four subsequent bytes
   PutPtr  -- Store a pointer into subsequent bytes

        BPtr    - Destination pointer
        WVal    - Word/Ptr

   This macro stores WVal into subsequent bytes. The pointer is advanced
   to the next free byte.

<*/
#define         PutWord(BPtr,WVal)                                      \
                {                                                       \
                   Take.Word = WVal;                                    \
                   PushWord (BPtr);                                     \
                }
#if PTR_SIZE == 2
#define         PutPtr(BPtr,WVal)                                       \
                {                                                       \
                   Take.Ptr = WVal;                                     \
                   *BPtr++ = Take.Byte [0];                             \
                   *BPtr++ = Take.Byte [1];                             \
                }
#else
#if PTR_SIZE == 4
#define         PutPtr(BPtr,WVal)                                       \
                {                                                       \
                   Take.Ptr = WVal;                                     \
                   *BPtr++ = Take.Byte [0];                             \
                   *BPtr++ = Take.Byte [1];                             \
                   *BPtr++ = Take.Byte [2];                             \
                   *BPtr++ = Take.Byte [3];                             \
                }
#else
#define         PutPtr(BPtr,WVal)                                       \
                {                                                       \
                   Take.Ptr = WVal;                                     \
                   *BPtr++ = Take.Byte [0];                             \
                   *BPtr++ = Take.Byte [1];                             \
                   *BPtr++ = Take.Byte [2];                             \
                   *BPtr++ = Take.Byte [3];                             \
                   *BPtr++ = Take.Byte [4];                             \
                   *BPtr++ = Take.Byte [5];                             \
                   *BPtr++ = Take.Byte [6];                             \
                   *BPtr++ = Take.Byte [7];                             \
                }
#endif
#endif
/*>

   GetWord -- Get a word from four subsequent bytes
   GetPtr  -- Get a pointer from subsequent bytes

        BPtr    - Byte pointer

   This macro assemblies word/pointer from subsequent bytes. The pointer
   is advanced to the next free byte.

<*/
#define         GetWord(BPtr)                                           \
                (                                                       \
                   PickUpWord (BPtr),                                   \
                   Take.Word                                            \
                )
#if PTR_SIZE == 2
#define         GetPtr(BPtr)                                            \
                (                                                       \
                   Take.Byte [0] = *BPtr++,                             \
                   Take.Byte [1] = *BPtr++,                             \
                   Take.Ptr                                             \
                )
#else
#if PTR_SIZE == 4
#define         GetPtr(BPtr)                                            \
                (                                                       \
                   Take.Byte [0] = *BPtr++,                             \
                   Take.Byte [1] = *BPtr++,                             \
                   Take.Byte [2] = *BPtr++,                             \
                   Take.Byte [3] = *BPtr++,                             \
                   Take.Ptr                                             \
                )
#else
#define         GetPtr(BPtr)                                            \
                (                                                       \
                   Take.Byte [0] = *BPtr++,                             \
                   Take.Byte [1] = *BPtr++,                             \
                   Take.Byte [2] = *BPtr++,                             \
                   Take.Byte [3] = *BPtr++,                             \
                   Take.Byte [4] = *BPtr++,                             \
                   Take.Byte [5] = *BPtr++,                             \
                   Take.Byte [6] = *BPtr++,                             \
                   Take.Byte [7] = *BPtr++,                             \
                   Take.Ptr                                             \
                )
#endif
#endif
/*>

   GetInt -- Assembling a signed integer from four subsequent bytes
   PutInt -- Storing a signed integer into four subsequent bytes

        Pointer - To the first byte

   The  Pointer is advanced so that after GetInt/PutInt it will point to
   four bytes far from initial position. The value is  stored/read  from
   Take.Int variable. The integer has portable format:

                Byte    Content
                 1      Sign (low significant bit) and bits 1..7
                 2      Bits 8..15
                 3      Bits 16..23
                 4      Bits 24..31

<*/
#define         GetInt(BPtr)                                            \
                {                                                       \
                   Take.Int  = 0x0000007FL & ((long) *BPtr++ >> 1);     \
                   Take.Int += 0x00007F80L & ((long) *BPtr++ << 7);     \
                   Take.Int += 0x007F8000L & ((long) *BPtr++ << 15);    \
                   Take.Int += 0x7F800000L & ((long) *BPtr++ << 23);    \
                   if (BPtr [-4] & 0x01) Take.Int = -Take.Int;          \
                }
#define         PutInt(BPtr)                                            \
                {                                                       \
                   if (0 > Take.Int)                                    \
                   {                                                    \
                      Take.Int = -Take.Int;                             \
                      *BPtr++ =						\
			 (byte) (((Take.Int & 0x7F) << 1) | 0x01);	\
                   }                                                    \
                   else                                                 \
                   {                                                    \
                      *BPtr++ = (byte) ((Take.Int & 0x7F) << 1);        \
                   }                                                    \
                   *BPtr++ = (byte) ((Take.Int >> 7 ) & 0xFF);          \
                   *BPtr++ = (byte) ((Take.Int >> 15) & 0xFF);          \
                   *BPtr++ = (byte) ((Take.Int >> 23) & 0xFF);          \
                }
/*>
        Statement header format.
        ------------------------
        This  format is used when a statement becomes active or matched.
        The purpose is to save some data  that  are  necessary  for  the
        statement processing. Any header contains a tag that defines its
        type. This tag occupies one byte. Therefore,  some  kludges  are
        used to avoid alignment of the C compiler. Noone can say that it
        is a good style of coding. It should be also mentioned that some
        compilers  always align data to the word boundary. Moreover some
        wonderful compilers have bugs  in  the  sizeof  function  for  a
        struct items. Function SizeOf is used to avoid these bugs.

<*/
#define         Hdr             Header_Tag
#define         Tag             Header_Tag
#define         SizeOf(Who)                                             \
                ((byte) (&Who.Header_Tag - (byte *) &Who + 1))

typedef union
{
#define         RANGE                   0
   struct               /* Pseudo statement                             */
   {
      byte              Tag;                    /* Header tag           */
   }  Range;

#define         SELECT                  1
   struct
   {
      byte *            Pat_ptr;                /* Pattern pointer      */
      byte              Tag;                    /* Header tag           */
   }  Select;

#define         NOT                     2
   struct
   {
      byte *            Pat_ptr;                /* Pattern pointer      */
      byte              Tag;                    /* Header tag           */
   }  Not;

#define         NOEMPTY                 3
   struct
   {
      word              Stub;                   /* Saved string offset  */
      byte              Tag;                    /* Header tag           */
   }  Noempty;

#define         REPEATER                4
   struct
   {
      byte *            Pat_ptr;                /* Pattern pointer      */
      word              Counter;                /* Repeatition count    */
      word              Wanted;                 /* Wanted count         */
      byte              Tag;                    /* Header tag           */
   }  Repeater;

#define         AS_MUCH                 5
   struct
   {
      word              Counter;                /* Repeatition count    */
      byte *            Pat_ptr;                /* Pattern pointer      */
      byte              Tag;                    /* Header tag           */
   }  As_much;

#define         AS_LITTLE               6
   struct
   {
      word              Counter;                /* Repeatition count    */
      byte *            Pat_ptr;                /* Pattern pointer      */
      byte              Tag;                    /* Header tag           */
   }  As_little;

#define         CALL                    7
   struct
   {
      byte *            Pat_ptr;                /* Pattern pointer      */
      byte              Tag;                    /* Header tag           */
   }  Call;

#define         ASSIGN_SUBSTR           8
   struct
   {
      word              Id;                     /* User provided id     */
      word              Offset;                 /* String offset        */
      byte              Tag;                    /* Header tag           */
   }  Assign;

#define         FENCED_ASSIGN           9
#define         SLEEPING_ASSIGN         10
   struct                                       /* Sleeping form        */
   {
      word              Id;                     /* User provided id     */
      byte              Tag;                    /* Header tag           */
   }  Sleeping_Assign;

#define         ELLIPSIS                11
   struct
   {
      byte *            Pat_ptr;                /* Pattern pointer      */
      byte              Tag;                    /* Header tag           */
   }  Ellipsis;
#define         TOTAL_STATEMENTS        12

                        /* This format is used for pattern translation  */
   struct
   {
      byte *            Name;                   /* Name string pointer  */
      word              Destination_offset;     /* Offset in pattern    */
   }  Label;

}       HeaderUnion;
/*>
        Stack manipulation

        PopStack        - poppes some data up from the stack
        FreeStack       - just the same, but without data moving
        PushStack       - pushes some data onto the stack
        GetHdr          - returns top stack byte
        SetStat         - takes header tag from stack and sets ThisStat
        StatTag         - makes header tag of the statement (ThisStat)

        In purpose of efficiency this code was made as macro.  This  has
        some disadvantages... !!WARNING!!

        1. Popping data must really be in the stack. No empty stack con-
           trol is made.
        2. Data may cross only one stack frame margin.
        3. The following variables must be visible (i.e.  declared)  and
           scratch: halfword DataLen, byte * DataPtr, byte * DataFrame
        4. Label NoMoreMemory must be visible too. You will be there  on
           malloc's failure.
<*/
#define         PopStack(Stack_ptr, Stack_used, Target, Size)           \
                {                                                       \
                   DataLen = (Size);                                    \
                   DataPtr = (byte *) (Target) + DataLen;		\
                   if (DataLen >= Stack_used)                           \
                   {                                                    \
                      DataLen = (halfword) (DataLen - Stack_used);      \
                      while (Stack_used--)                              \
                         *--DataPtr = *--Stack_ptr;                     \
                      PopPtr (Stack_ptr);				\
                      Stack_ptr = Take.Ptr;				\
                      Stack_used = MAX_STACK_FRAME;			\
		   }							\
                   Stack_used = (halfword) (Stack_used - DataLen);      \
                   while (DataLen--)                                    \
                      *--DataPtr = *--Stack_ptr;                        \
                }

#define         FreeStack(Stack_ptr, Stack_used, Size)                  \
                {                                                       \
                   DataLen = (Size);                                    \
                   if (DataLen >= Stack_used)                           \
                   {                                                    \
                      DataLen = (halfword) (DataLen - Stack_used);      \
                      Stack_ptr -= Stack_used;                          \
                      PopPtr (Stack_ptr);				\
                      Stack_ptr = Take.Ptr;				\
                      Stack_used = MAX_STACK_FRAME;			\
                   }							\
                   Stack_used = (halfword) (Stack_used - DataLen);      \
                   Stack_ptr -= DataLen;                                \
                }

#define         PushStack(Stack_ptr, Stack_used, Target, Size)          \
                {                                                       \
                   DataLen = (Size);                                    \
                   DataPtr = (byte *) (Target);                         \
                   if (DataLen > MAX_STACK_FRAME - Stack_used)          \
                   {                                                    \
                      Stack_used =				        \
		         (halfword) (MAX_STACK_FRAME - Stack_used);     \
                      DataLen = (halfword) (DataLen - Stack_used);      \
                      while (Stack_used--)                              \
                         *Stack_ptr++ = *DataPtr++;                     \
                      if (GetPtr (Stack_ptr) != 0)                      \
                      {                                                 \
                         Stack_ptr = Take.Ptr;                          \
                      }                                                 \
                      else                                              \
                      {                                                 \
                         if (  NULL                                     \
                            == (  DataFrame                             \
                                = (byte *) GC_malloc (STACK_FRAME_SIZE)    \
                            )  )  goto NoMoreMemory;                    \
                         PushPtr (Stack_ptr, DataFrame + PTR_SIZE);     \
                         PutPtr (DataFrame, Stack_ptr);			\
                         Stack_ptr = DataFrame;                         \
                         DataFrame += MAX_STACK_FRAME;                  \
                         PutPtr (DataFrame, 0);                         \
                      }                                                 \
                      Stack_used = 0;                                   \
                   }                                                    \
                   Stack_used = (halfword) (Stack_used + DataLen);      \
                   while (DataLen--)                                    \
                      *Stack_ptr++ = *DataPtr++;                        \
                }

#define         GetHdr(Pointer)         (((byte *) Pointer) [-1])

#define         SetStat(Pointer)                                        \
                (ThisStat = (HEADER_SLICE_MASK & GetHdr (Pointer)) >>1)

#define         StatTag                                                 \
                ((byte) (HEADER_TEST_MASK | (ThisStat<<1)))

/*>     Common Stack processing

        SaveStack       - saves current stack state (reentrant call)
        RestoreStack    - restores stack state
        StackData       - defines data to save stack state
        Done            - an equivalent to ordinal return
        ErrorDone       - set MatchError and Done
<*/
#define         SaveStack                                               \
                {                                                       \
                   First_ASS_ptr  = ASS_ptr;                            \
                   First_MSS_ptr  = MSS_ptr;                            \
                   First_ASS_used = ASS_used;                           \
                   First_MSS_used = MSS_used;                           \
                }

#define         RestoreStack                                            \
                {                                                       \
                   ASS_ptr  = Old_ASS_ptr;                              \
                   MSS_ptr  = Old_MSS_ptr;                              \
                   ASS_used = Old_ASS_used;                             \
                   MSS_used = Old_MSS_used;                             \
                   First_ASS_ptr  = Old_1st_ASS_ptr;                    \
                   First_MSS_ptr  = Old_1st_MSS_ptr;                    \
                   First_ASS_used = Old_1st_ASS_used;                   \
                   First_MSS_used = Old_1st_MSS_used;                   \
                }

#define         StackData                                               \
                   byte *       Old_ASS_ptr  = ASS_ptr;                 \
                   byte *       Old_MSS_ptr  = MSS_ptr;                 \
                   halfword     Old_ASS_used = ASS_used;                \
                   halfword     Old_MSS_used = MSS_used;                \
                   byte *       Old_1st_ASS_ptr  = First_ASS_ptr;       \
                   byte *       Old_1st_MSS_ptr  = First_MSS_ptr;       \
                   halfword     Old_1st_ASS_used = First_ASS_used;      \
                   halfword     Old_1st_MSS_used = First_MSS_used;

#define         Done(Expression)                                        \
                {                                                       \
                   RestoreStack;                                        \
                   return (Expression);                                 \
                }

#define         ErrorDone(Expression)                                   \
                {                                                       \
                   RestoreStack;                                        \
                   return ( MatchError = (Expression) );                \
                }

/*>     Active Statement Stack processing

        EmptyASS        - returns TRUE if ASS is empty
        PopASS          - takes the statement (ThisStat) header from ASS
        PushASS         - pushes the statement header onto ASS
        PopLabel        - pop a label descriptor from ASS (Patran only)
        PushLabel       - push a label descriptor onto ASS
<*/
#define         EmptyASS                                                \
                (ASS_ptr == First_ASS_ptr && ASS_used == First_ASS_used)

#define         PopASS                                                  \
                PopStack                                                \
                (  ASS_ptr, ASS_used,                                   \
                   (byte *) &Header,                                    \
                   StatSize [ThisStat]                                  \
                )

#define         PushASS                                                 \
                PushStack                                               \
                (                                                       \
                   ASS_ptr, ASS_used,                                   \
                   (byte *) &Header,                                    \
                   StatSize [ThisStat]                                  \
                )

#define         PopLabel                                                \
                PopStack                                                \
                (                                                       \
                   ASS_ptr, ASS_used,                                   \
                   (byte *) &Header.Label,                              \
                   sizeof (Header.Label)                                \
                )

#define         PushLabel                                               \
                PushStack                                               \
                (                                                       \
                   ASS_ptr, ASS_used,                                   \
                   (byte *) &Header.Label,                              \
                   sizeof (Header.Label)                                \
                )

/*>     Matched Statement Stack processing

        EmptyMSS        - returns TRUE if MSS is empty
        PopMSS          - takes the statement (ThisStat) header from MSS
        PushMSS         - pushes the statement header onto MSS
        StubOnMSS       - returns TRUE if a stub lies on the top of MSS
        KillStub        - just removes a stub from the MSS
        CheckString     - order Str_ptr if Length and Pointer are incorrect
        SaveStr         - saves string as a stub on the MSS
        RestoreStr      - restores string pointer from the stub
        MakeStub        - makes a stub value
<*/
#define         EmptyMSS                                                \
                (MSS_ptr == First_MSS_ptr && MSS_used == First_MSS_used)

#define         PopMSS                                                  \
                PopStack                                                \
                (                                                       \
                   MSS_ptr, MSS_used,                                   \
                   (byte *) &Header,                                    \
                   StatSize [ThisStat]                                  \
                )

#define         PushMSS                                                 \
                PushStack                                               \
                (                                                       \
                   MSS_ptr, MSS_used,                                   \
                   (byte *) &Header,                                    \
                   StatSize [ThisStat]                                  \
                )

#define         StubOnMSS                                               \
                (0 == (GetHdr (MSS_ptr) & 0x01))

#define         KillStub        FreeStack (MSS_ptr, MSS_used, WORD_SIZE)

#define         MakeStub						\
		(((word) (Str_ptr - Str_Beg) + Str_offset) << 1)

#define         CheckString                                             \
                {                                                       \
                   if (Length < 0) Str_End = Str_Beg;                   \
                   if (*Pointer < 0)                                    \
                      Str_ptr = Str_Beg;                                \
                   else                                                 \
                      if (*Pointer > Length)                            \
                         Str_ptr = Str_End;                             \
                }
#define         SaveStr                                                 \
                {                                                       \
		   Temp = MakeStub;					\
                   DataLen = WORD_SIZE;					\
                   if (DataLen > MAX_STACK_FRAME - MSS_used)		\
                   {                                                    \
                      MSS_used =				        \
		         (halfword) (MAX_STACK_FRAME - MSS_used);	\
                      DataLen = (halfword) (DataLen - MSS_used);	\
                      switch (MSS_used)					\
		      {							\
		         case 4 : MSS_ptr [3] =(byte)(0xFF &(Temp    ));\
		         case 3 : MSS_ptr [2] =(byte)(0xFF &(Temp>> 8));\
		         case 2 : MSS_ptr [1] =(byte)(0xFF &(Temp>>16));\
		         case 1 : MSS_ptr [0] =(byte)(0xFF &(Temp>>24));\
		      }							\
		      MSS_ptr += MSS_used;				\
                      if (GetPtr (MSS_ptr) != 0)			\
                      {                                                 \
                         MSS_ptr = Take.Ptr;				\
                      }                                                 \
                      else                                              \
                      {                                                 \
                         if (  NULL                                     \
                            == (  DataFrame                             \
                                = (byte *) GC_malloc (STACK_FRAME_SIZE)    \
                            )  )  goto NoMoreMemory;                    \
                         PushPtr (MSS_ptr, DataFrame + PTR_SIZE);	\
                         PutPtr (DataFrame, MSS_ptr);			\
                         MSS_ptr = DataFrame;				\
                         DataFrame += MAX_STACK_FRAME;                  \
                         PutPtr (DataFrame, 0);                         \
                      }                                                 \
                      MSS_used = 0;					\
                   }                                                    \
                   MSS_used = (halfword) (MSS_used + DataLen);		\
                   MSS_ptr += DataLen;					\
                   switch (DataLen)					\
		   {							\
		      case 4 : MSS_ptr [-4] = (byte)(0xFF &(Temp>>24));	\
		      case 3 : MSS_ptr [-3] = (byte)(0xFF &(Temp>>16));	\
		      case 2 : MSS_ptr [-2] = (byte)(0xFF &(Temp>> 8));	\
		      case 1 : MSS_ptr [-1] = (byte)(0xFF &(Temp    ));	\
                   }							\
                }
#define         RestoreStr                                              \
                {                                                       \
		   Temp = 0;						\
		   DataLen = WORD_SIZE;					\
                   if (DataLen >= MSS_used)				\
                   {							\
		      DataLen = (halfword) (DataLen - MSS_used);	\
		      switch (MSS_used)					\
		      {							\
		         case 4 : Temp += (word) MSS_ptr [-4] << 24;	\
		         case 3 : Temp += (word) MSS_ptr [-3] << 16;	\
		         case 2 : Temp += (word) MSS_ptr [-2] << 8;	\
		         case 1 : Temp += (word) MSS_ptr [-1];		\
		      }							\
                      MSS_ptr -= MSS_used;	                        \
                      PopPtr (MSS_ptr);					\
                      MSS_ptr = Take.Ptr;				\
                      MSS_used = MAX_STACK_FRAME;			\
                   }							\
                   MSS_used = (halfword) (MSS_used - DataLen);		\
                   MSS_ptr -= DataLen;					\
		   switch (DataLen)					\
		   {							\
		      case 4 : Temp += (word) MSS_ptr [3];		\
		      case 3 : Temp += (word) MSS_ptr [2] << 8;		\
		      case 2 : Temp += (word) MSS_ptr [1] << 16;	\
		      case 1 : Temp += (word) MSS_ptr [0] << 24;	\
		   }							\
		   if (Str_offset > (Temp >>= 1))			\
                   {                                                    \
		      int   Length;					\
									\
                      while (Temp < Str_offset)                         \
                      {                                                 \
                         if (  GetPreviousLine == NULL                  \
                            || ! (* GetPreviousLine)                    \
                                 ( (const char **) &Str_Beg, &Length	\
                            )    )  goto CantReturn;                    \
                         Str_offset -= Length + 1;			\
                      }                                                 \
                      Str_End = Str_Beg + Length;			\
                   }                                                    \
                   Str_ptr = Str_Beg + (natural) (Temp - Str_offset);   \
                }

/*>     Source string processing

        DoCircumflex    - byte conversion by a circumflex
        EndOfString     - returns TRUE if end of the string seen
        Error           - error logging (fatal)
        IsUserAlpha     - is character is user defined alphabeticals
        KeyWord         - returns TRUE if specified keyword is matched
        PassSpace       - skip spaces and tabs
        PassId          - skip identifier

<*/
#define         DoCircumflex(Code)                                      \
                ((Code) ^= 0100)

#define         EndOfString             (Str_ptr >= Str_End)

#define         IsUserAlpha(Char)                                       \
                (  isalpha (Code = (Char))                              \
                || (  UserAlphas != NULL                                \
                   && strchr (UserAlphas, (Code)) != NULL               \
                )  )

#define         KeyWord(Text)                                           \
                (  Lex_ptr						\
		!= (  Str_ptr						\
		   =  MatchKeyWord (Lex_ptr, Str_End, (Text))		\
                )  )

#define         PassSpace                                               \
                {                                                       \
                   while (  !EndOfString                                \
                         && ( (Code = *Str_ptr) == SPACE                \
                            || Code == TAB                              \
                         )  )  Str_ptr++;                               \
                }

#define         PassId                                                  \
                {                                                       \
                   while (  !EndOfString                                \
                         && (  IsUserAlpha (*Str_ptr)                   \
                            || isdigit (Code)                           \
                         )  )  Str_ptr++;                               \
                   Lex_End = Str_ptr;                                   \
                   PassSpace;                                           \
                }

#define         Error(Code, Pointer)                                    \
                {                                                       \
                   *Error_code = Code;                                  \
                   *Error_ptr  = Pointer;                               \
                   return (0);                                          \
                }

/*>     Brace processing

        PopBrace        - pop a brace (resulted in 0 - (), 1 - [] or 2)

        !!WARNING!! Variable Zero must contain  0  and  be  visible  for
                    PushXB.

<*/
#define         PopBrace(Result)                                        \
                {                                                       \
                   SAY (4, "\n           ^ MSS_mask = %x", MSS_mask);   \
                   if (MSS_mask == 0)                                   \
                   {                                                    \
                      Result = 2;                                       \
                   }                                                    \
                   else                                                 \
                   {                                                    \
                      Result =						\
		         (byte)	(0 != (*(MSS_ptr - 1) & MSS_mask));	\
                      if ((MSS_mask >>= 1) == 0)                        \
                      {                                                 \
                         FreeStack (MSS_ptr, MSS_used, 1);              \
                         if (!EmptyMSS) MSS_mask = LAST_MSS_MASK;       \
                }  }  }

#ifdef DEBUG
#define         SAY(level, format, data)                                \
                if (DEBUG > level) fprintf (stderr, (format), (data))
#else
#define         SAY(level, format, data)
#endif

/*>     Statement codes to detect indefinite loops (see Second_pass)    <*/

#define         EMPTY_ATOM        (     0     ) /* matches null string  */
#define         NON_EMPTY_ATOM    (     1     ) /* does not match null  */
#define         NON_SELECT_MASK   (     2     ) /* non-select if set    */
#define         RETURN_EMPTY      ( 2+( 0<<2 )) /* stat. can match null */
#define         RETURN_NON_EMPTY  ( 2+( 1<<2 )) /* stat. can't do it    */
#define         SON_WILL_SAY      ( 2+( 2<<2 )) /* depending on subpat  */
#define         INDEFINITE        ( 2+( 3<<2 )) /* $ or * repeater      */
#define         SELECT_SHIFT      (     4     ) /* shift to empty flag  */
#define         SELECT_MASK       (~( 1<<4 )-1) /* mask for empty flag  */
#define         SELECT_BRANCH     (    3<<4   ) /* select statement     */
#define         ERROR_STATUS      (    2<<5   ) /* error flag           */

/*>     Macros ...

        PushStatement   - push a statement code
        PopStatement    - pop it
        FindForSelect   - pop until select and test empty branch existing
        CloseBranch     - accumulate branch status
<*/

#define         PushStatement                                           \
                PushStack (MSS_ptr, MSS_used, &EmptyStatus, 1)

#define         PopStatement                                            \
                PopStack (MSS_ptr, MSS_used, &EmptyStatus, 1);

#define         FindForSelect(Status)                                   \
                if (  (  EmptyStatus = CloseStatement (Status)          \
                      )  == ERROR_STATUS                                \
                   )  goto IndefiniteLoop;

#define         CloseBranch(Status)                                     \
                FindForSelect(Status);                                  \
                *(MSS_ptr-1) &= (  ( *(MSS_ptr-1) << SELECT_SHIFT )     \
                                |  SELECT_MASK                          \
                                );

/*^v^v^v^v^v^v^v^ STATIC AND COMMON VARIABLES ^v^v^v^v^v^v^v^v^v^v^v^v^*/

static byte     PowerOf2 [8]    = { 1, 2, 4, 8, 16, 32, 64, 128 };
static byte     Preallocated_ASS_frame [STACK_FRAME_SIZE];
static byte     Preallocated_MSS_frame [STACK_FRAME_SIZE];
static byte *   ASS_ptr         = &Preallocated_ASS_frame [PTR_SIZE];
static byte *   First_ASS_ptr   = &Preallocated_ASS_frame [PTR_SIZE];
static byte *   First_MSS_ptr   = &Preallocated_MSS_frame [PTR_SIZE];
static byte *   MSS_ptr         = &Preallocated_MSS_frame [PTR_SIZE];
static halfword First_ASS_used  = 1;	/* to prevent deallocation of	*/
static halfword First_MSS_used  = 1;	/* the preallocated frames	*/
static halfword ASS_used        = 1;    /* the ASS frame used bytes     */
static halfword MSS_used        = 1;    /* the MSS frame used bytes     */
static halfword MSS_mask;               /* the bit mask of MSS top byte */
static byte *   DataPtr;
static byte *   DataFrame;
static halfword DataLen;

#define         Shared
#ifndef NON_ANSI
int           (* GetNextLine)        (const char ** Line, int * Length) = 0;
int           (* GetPreviousLine)    (const char ** Line, int * Length) = 0;
char *	      (* AllocatePattern)    (unsigned long Size)		= 0;
char *        (* GetExternalPattern) (char *        Name, int   Length) = 0;
// unsigned long (* GetVariableId)      (char *        Name, int   Length)	= 0;
/*void          (* AssignVariable)     (  unsigned long VariableId,
                                        unsigned long Offset,
                                        unsigned long Length
                                     )					= 0;*/
// void          (* DeAssignVariable)   (unsigned long VariableId)		= 0;
#else
int           (* GetNextLine)        () = 0;
int           (* GetPreviousLine)    () = 0;
char *	      (* AllocatePattern)    () = 0;
char *        (* GetExternalPattern) () = 0;
unsigned long (* GetVariableId)      () = 0;
void          (* AssignVariable)     () = 0;
void          (* DeAssignVariable)   () = 0;
#endif /*> NON_ANSI <*/

char *  UserAlphas      = NULL; /* User defined alphabeticals   */
int     UserBase        = 10;   /* Digits are 0..9              */
int     MatchError      = 0;    /* Status code                  */

static struct
   {
      byte      Code;
      char *    Name;
   }
   KeyList [KEY_LIST_SIZE] =
   {
      {  ANY               , "ANY"        },
      {  ANY               , "%"          },
      {  ARB               , "ARB"        },
      {  ARB               , "..."        },
      {  ARB               , ".."         },
      {  BIG_REPEATER      , "$"          },
      {  BLANK             , "BLANK"      },
      {  BLANK             , "+"          },
      {  CHARACTER         , "CHARACTER"  },
      {  CHARACTER         , "C"          },
      {  DIGIT             , "DIGIT"      },
      {  DIGIT             , "#"          },
      {  DO_NOT_RETURN     , "FENCE"      },
      {  DO_NOT_RETURN     , ":"          },
      {  END_OF_KEYWORD    , "BREAK"      },
      {  END_OF_KEYWORD    , "_"          },
      {  END_OF_STRING     , "END"        },
      {  END_OF_STRING     , "."          },
      {  FAILURE           , "FAILURE"    },
      {  FAILURE           , "F"          },
      {  INVERSE           , "NOT"        },
      {  INVERSE           , "^"          },
      {  LETTER            , "LETTER"     },
      {  LOWER_CASE_LETTER , "LCL"        },
      {  LETTER            , "L"          },
      {  LITTLE_REPEATER   , "*"          },
      {  LOWER_CASE_LETTER , "W"          },
      {  NEW_LINE          , "NL"         },
      {  NEW_LINE          , "/"          },
      {  OR                , "OR"         },
      {  OR                , "!"          },
      {  OR                , "|"          },
      {  QUERY             , "NOEMPTY"    },
      {  QUERY             , "?"          },
      {  SUCCESS           , "SUCCESS"    },
      {  SUCCESS           , "S"          },
      {  UPPER_CASE_LETTER , "UCL"        },
      {  UPPER_CASE_LETTER , "U"          },
   };

static byte     StatSize [TOTAL_STATEMENTS];    /* Sizes of headers     */

static union
{
   long         Int;
   word         Word;
   byte *       Ptr;
   byte         Byte [PTR_SIZE];
}               Take;

HeaderUnion     Header;

/*^v^v^v^v^v^v^v^v^v^v^v^ PSEUDO STATEMENTS ^v^v^v^v^v^v^v^v^v^v^v^v^v^*/

#define         SkipOneItem_and_Success                                 \
                { ReturnTo = 0 ; goto SkipOneItem; }

#define         SkipOneItem_and_IsAlternative                           \
                { ReturnTo = 1 ; goto SkipOneItem; }

#define         SkipOneItem_and_EndOfSelect                             \
                { ReturnTo = 2 ; goto SkipOneItem; }

#define         SkipOneItem_and_DoLittleRepeater                        \
                ReturnTo = 3 ; goto SkipOneItem; DoLittleRepeater :

#define         CleanMSS_and_Success                                    \
                ReturnTo = 0 ; goto CleanMSS;

#define         CleanMSS_and_Failure                                            \
                ReturnTo = 1 | 0x10 | 0x20; goto CleanMSS;

#define         CleanMSS_and_Go                                         \
                ReturnTo = 2 | 0x10; goto CleanMSS;

#define         CleanMSS_and_ContinueLittleRepeater                     \
                ReturnTo = 3 | 0x10; goto CleanMSS;                     \
                ContinueLittleRepeater :

/*^v^v^v^v^v^v^v^v^v^v^v^v^v^v^ FUNCTIONS ^v^v^v^v^v^v^v^v^v^v^v^v^v^v^v

   External functions

   freedm   -- releases dynamic memory occupied by working stacks
   match    -- matches string with a pattern
   patran   -- pattern translator
   patmaker -- pattern constructor


   freedm -- releases dynamic memory occupied by ASS and MSS

   After execution ASS and MSS pointers will  be  restored  to  previous
   values  (First_ASS_ptr,  First_ASS_used, ...). All stack frames allo-
   cated later will be released.

                              Will be released >>>>>>>>>>>>>>>>>>>>>>>
                ._________.      ._________.      ._________.
                |         |      |         |      |         |
            --->|  frame  |----->|  frame  |----->|  frame  |----> nil
                |_________|      |_________|      |_________|
                       |              |
       First_xSS_ptr __|    xSS_ptr __|                     x = A or M

<*/
#ifndef NON_ANSI
void freedm (void)
#else
void freedm ()
#endif                  /*> NON_ANSI <*/
{
  register byte *      ThisBlock;
  register byte *      NextBlock;
  int                  StackNo;

  for (StackNo = 0; StackNo < 2; StackNo++)
  {
    ThisBlock =
        (  (  StackNo
        ?  First_ASS_ptr - First_ASS_used
      :  First_MSS_ptr - First_MSS_used
           )
        +  MAX_STACK_FRAME
        );
    NextBlock = GetPtr (ThisBlock);
    PushPtr (ThisBlock, 0);
    while (NextBlock != 0)
    {
      ThisBlock = NextBlock + MAX_STACK_FRAME;
      NextBlock = GetPtr (ThisBlock);
      ThisBlock -= STACK_FRAME_SIZE;
      GC_free (ThisBlock);
    }
  }
  ASS_ptr  = First_ASS_ptr;
  MSS_ptr  = First_MSS_ptr;
  ASS_used = First_ASS_used;
  MSS_used = First_MSS_used;
}       /*>  freedm  <*/
/*>

   match -- match a pattern

        Length  - of a character string to be matched
        String  - points to the string (maybe not NUL terminated)
        Pointer - number of 1st character to be matched (0..Length)
        Pattern - points to a pattern (in internal format)

   Returns :
                     MATCH_SUCCESS - Success
                     MATCH_FAILURE - Failure of matching
        MATCH_WRONG_PATTERN_FORMAT - Wrong pattern format
             MATCH_TOO_LONG_STRING - String length exceeds 1Gbyte
               MATCH_CANNOT_RETURN - Cannot return to previously left string
           MATCH_NO_DYNAMIC_MEMORY - Dynamic memory overflow

   After successful completion String [*Pointer] will be the 1st charac-
   ter  following matched substring. Note that after successful matching
   the new line atom (/) Pointer will be still counted as an offset from
   beginning  of  the String. The offset includes line ends: +1 for each
   end.

<*/

#ifndef NON_ANSI

int match
(
   int          Length,
   const char * String,
   int        * Pointer,
   const char * Pattern
)
#else

int match (Length, String, Pointer, Pattern)

   int          Length;
   char       * String;
   int        * Pointer;
   char       * Pattern

#endif  /*> NON_ANSI <*/
{
   register byte *      Str_ptr = (byte *) &String [*Pointer];
   register byte *      Pat_ptr = (byte *) Pattern;
   register byte        Code;
   byte          *      Str_End = (byte *) &String [Length];
   byte          *      Str_Beg = (byte *) String;
   byte                 ReturnTo;
   word                 Str_offset = 0;
   word                 Temp;
   int			ThisStat;
   HeaderUnion          Header;

   StackData;

#ifdef DEBUG
   assert (SizeOf (Header.Select) == sizeof (Header.Select.Pat_ptr) + 1);
#endif /*> DEBUG <*/

   /*> Something to initialize <*/

   if (StatSize [RANGE] == 0)
   {
      StatSize [RANGE]           = SizeOf (Header.Range);
      StatSize [SELECT]          = SizeOf (Header.Select);
      StatSize [NOT]             = SizeOf (Header.Not);
      StatSize [NOEMPTY]         = SizeOf (Header.Noempty);
      StatSize [REPEATER]        = SizeOf (Header.Repeater);
      StatSize [AS_MUCH]         = SizeOf (Header.As_much);
      StatSize [AS_LITTLE]       = SizeOf (Header.As_little);
      StatSize [CALL]            = SizeOf (Header.Call);
      StatSize [ASSIGN_SUBSTR]   = SizeOf (Header.Assign);
      StatSize [FENCED_ASSIGN]   = SizeOf (Header.Sleeping_Assign);
      StatSize [SLEEPING_ASSIGN] = SizeOf (Header.Sleeping_Assign);
      StatSize [ELLIPSIS]        = SizeOf (Header.Ellipsis);
   }
#ifdef DEBUG
   assert (SizeOf (Header.Select) == sizeof (Header.Select.Pat_ptr) + 1);
   assert (MAX_STACK_FRAME > StatSize [RANGE]);
   assert (MAX_STACK_FRAME > StatSize [SELECT]);
   assert (MAX_STACK_FRAME > StatSize [NOT]);
   assert (MAX_STACK_FRAME > StatSize [NOEMPTY]);
   assert (MAX_STACK_FRAME > StatSize [REPEATER]);
   assert (MAX_STACK_FRAME > StatSize [AS_MUCH]);
   assert (MAX_STACK_FRAME > StatSize [AS_LITTLE]);
   assert (MAX_STACK_FRAME > StatSize [CALL]);
   assert (MAX_STACK_FRAME > StatSize [ASSIGN_SUBSTR]);
   assert (MAX_STACK_FRAME > StatSize [FENCED_ASSIGN]);
   assert (MAX_STACK_FRAME > StatSize [SLEEPING_ASSIGN]);
   assert (MAX_STACK_FRAME > StatSize [ELLIPSIS]);
   assert (MAX_STACK_FRAME > WORD_SIZE);
#endif /*> DEBUG <*/
   CheckString;
   SaveStack;
   ThisStat = RANGE;
   Header.Range.Hdr = StatTag;
   PushASS;
   PushMSS;
   ThisStat = SELECT;
   Header.Select.Hdr = StatTag;
   Header.Select.Pat_ptr = Pat_ptr;
   SAY (4, "\nMatch      Machine word of %d bytes", MachineWord);

Activate :      /*>     Here a statement indicated by the  ThisStat  be-
                        comes  active.  Its  header is pushed to the ASS
                        and current string pointer is saved in the MSS.
                                                                        <*/
   PushASS;

Go :            /*>     Save string and take the next item of the
                        pattern                                         <*/
   SaveStr;

Process :

   Code = *Pat_ptr++;
   SAY (4, "\nMatch      PROCESS %o (Oct)", Code);

   if (Code <= LAST_ASCII_CHARACTER)
   {/*>         Match a 7-bit ASCII text constant       <*/
      SAY (2, "\nMatch      Do literal %c", Code);
      while (Code != SUCCESS)
      {
         if (EndOfString || Code != *Str_ptr++) goto Failure;
         Code = *Pat_ptr++;
         SAY (2, "%c", Code);
         if (Code > LAST_ASCII_CHARACTER)
         {
            --Pat_ptr;
            goto Success;
      }  }
      SAY (2, "\nMatch      GENERAL SUCCESS", Code);
      goto GeneralSuccess;              /* Atom 'success' seen          */
   }
   if ((natural) Code <= MAX_FINITE_REPEATER)
   {
      Header.Repeater.Wanted = Code - LAST_ASCII_CHARACTER - 1;

MakeRepeater :

      SAY (2, "\nMatch      Do FINITE REPEATER for %ld times", Header.Repeater.Wanted);
      if (Header.Repeater.Wanted)
      {
         ThisStat = REPEATER;
         Header.Repeater.Counter = 0;
         Header.Repeater.Pat_ptr = Pat_ptr;
         Header.Repeater.Hdr = StatTag;
         PushASS;
         goto Process;
      }
      else
      {
         SkipOneItem_and_Success;
   }  }

   if ((natural) Code < FIRST_STATEMENT_CODE)
   {
      SAY (2, "\nMatch      Do atom %o (Oct) ", Code);
      switch (Code)
      {
         case END_OF_STRING :
            SAY (2, " - End of string", Code);
            if (EndOfString) goto Success;
            goto Failure;
         case END_OF_KEYWORD :
            SAY (2, " - End of keyword", Code);
            if (EndOfString) goto Success;
            Code = *Str_ptr++;
            if (Code == SPACE || Code == TAB)
            {
               PassSpace;
               goto Success;
            }
            if (--Str_ptr == Str_Beg) goto Success;
            if (isalnum (Code) && isalnum (*(Str_ptr-1))) goto Failure;
         case NOOP :
            goto Success;
         case NEW_LINE :
            {
               char *   NewLine;
               int      Length;

               SAY (2, " - New Line", Code);
               if (  GetNextLine == NULL
                  || ! (* GetNextLine) ((const char **) &NewLine, &Length)
                  )  goto Failure;
               Temp = Str_End - Str_Beg + 1;
               if (Str_offset > MAX_STRING_LENGTH - Temp)
               {
                  ErrorDone (MATCH_TOO_LARGE_STRING);
               }
               Str_offset += Temp;
               Str_Beg = (byte *) NewLine;
               Str_End = Str_Beg + Length;
               Str_ptr = Str_Beg;
               goto Success;
            }
         case FAILURE :
            ErrorDone (MATCH_FAILURE);
      }
      SAY (2, " - ", Code);
      if (EndOfString) goto Failure;
      switch (Code)
      {
         case ANY :
            SAY (2, "ANY", Code);
            Str_ptr++;
            goto Success;
         case BLANK :
            SAY (2, "BLANK", Code);
            Code = *Str_ptr++;
            if (Code != SPACE && Code != TAB) goto Failure;
            PassSpace;
            goto Success;
         case DIGIT :
            SAY (2, "DIGIT", Code);
            Code = *Str_ptr++;
            if (UserBase <= 10)
            {  /* Bases 1..10                                   */
               if (  Code >= '0'
                  && Code <  '0' + UserBase
                  )  goto Success;
            }
            else
            {  /* Bases 11..                                    */
               if (  (Code >= '0' && Code <= '9')
                  || (  (  Code = (byte) tolower (Code)
                        )  >= 'a'
                     && Code <= 'a' + UserBase - 11
                  )  )  goto Success;
            }
            goto Failure;
         case UPPER_CASE_LETTER :
            SAY (2, "UPPER CASE LETTER", Code);
            Code = *Str_ptr++;
            if (isupper (Code)) goto Success;
            goto Failure;
         case LOWER_CASE_LETTER :
            SAY (2, "LOWER CASE LETTER", Code);
            Code = *Str_ptr++;
            if (islower (Code)) goto Success;
            goto Failure;
         case LETTER :
            SAY (2, "LETTER", Code);
            Code = *Str_ptr++;
            if (isalpha (Code)) goto Success;
            goto Failure;
         case CHARACTER :
            SAY (2, "CHARACTER", Code);
            Code = *Str_ptr++;
            if (isalnum (Code)) goto Success;
            goto Failure;
         case ANY_OF :
         {
            int Index, Size;

            SAY (2, "ANY OF", Code);
            if (  (  0
                  >  (  Index
                     =  ((int) (Code = *Str_ptr++) >> 3) - *Pat_ptr++
                  )  )
               || Index >= (Size = *Pat_ptr++)
               || 0 == (Pat_ptr [Index] & PowerOf2 [Code & 0x07])
               )  goto Failure;
            Pat_ptr += Size;
            goto Success;
         }
         case HOLERITH :
            {
               register int     Length;

               SAY (2, "HOLERITH's literal", Code);
               Length = 1 + (int) *Pat_ptr++;
               if ((Str_ptr + Length) > Str_End) goto Failure;
               while (Length--)
                  if (*Str_ptr++ != *Pat_ptr++) goto Failure;
               goto Success;
            }
         case CASE_DEAF_HOLERITH :
            {
               register int     Length;

               SAY (2, "CASE DEAF HOLERITH's literal", Code);
               Length = 1 + (int) *Pat_ptr++;
               if ((Str_ptr + Length) > Str_End) goto Failure;
               while (Length--)
                  if (tolower (*Str_ptr++) != *Pat_ptr++) goto Failure;
               goto Success;
            }
         default :                      /* Any other leads to failure   */
            ErrorDone (MATCH_WRONG_PATTERN_FORMAT);
   }  }
   switch (Code)
   {
      case EXTERNAL_PATTERN :
         SAY (2, "\nMatch      Do absolute CALL", Code);
         ThisStat = CALL;
         Header.Call.Hdr = StatTag;
         Header.Call.Pat_ptr = Pat_ptr + PTR_SIZE;
         Pat_ptr = GetPtr (Pat_ptr);
         SAY (8, "\nMatch      New Pat_ptr = %lx (Hex)", Pat_ptr);
         SAY (2, "\n           Now, emulate a left brace occurence", Code);
         SaveStr;
         PushASS;
      case LEFT_BRACE :
         SAY (2, "\nMatch      Do left brace", Code);
         ThisStat = SELECT;
         Header.Select.Hdr = StatTag;
         Header.Select.Pat_ptr = Pat_ptr;
         goto Activate;
      case RIGHT_BRACE :
         SetStat (ASS_ptr);
         if (ThisStat != SELECT)
         {
            --Pat_ptr;
            goto Success;
         }                      /* It's true select statement           */
                                /* We'll member current pattern pointer */
         PopASS;                /* to optimize next alternative search  */
         Header.Select.Pat_ptr = Pat_ptr - 1;
         SAY (2, "\nMatch      Do right brace", Code);
         goto Sleep;
      case OR :
         SetStat (ASS_ptr);
         if (ThisStat != SELECT)
         {
            --Pat_ptr;
            goto Success;
         }                      /* It's true select statement           */
                                /* We'll member current pattern pointer */
         PopASS;                /* to optimize next alternative search  */
         Header.Select.Pat_ptr = Pat_ptr - 1;
         SAY (2, "\nMatch      Do |", Code);

EndOfSelect :                   /* Find for end of the select statement */

         switch (*Pat_ptr)
         {
            case RIGHT_BRACE :  /* Ok, now we reach the end of select   */
               Pat_ptr++;       /* Skip the right brace and continue    */
               goto Sleep;
            case OR :
               Pat_ptr++;       /* Skip '|' and the next item too       */
            default :
               SkipOneItem_and_EndOfSelect;
         }
      case ARB :
         SAY (2, "\nMatch      Do ... - arbitrary number of chars", Code);
         SaveStr;
         ThisStat = ELLIPSIS;
         Header.Ellipsis.Hdr     = StatTag;
         Header.Ellipsis.Pat_ptr = Pat_ptr;
         goto Sleep;
      case FINITE_REPEATER :
         GetInt (Pat_ptr);
         Header.Repeater.Wanted = Take.Int;
         goto MakeRepeater;
      case GO_TO :
         ThisStat = CALL;
         Header.Call.Hdr = StatTag;
         GetInt (Pat_ptr);
         SAY (2, "\nMatch      Do relative CALL at %ld", Take.Int);
         Header.Call.Pat_ptr = Pat_ptr;
         Pat_ptr += Take.Int;
         goto Activate;
      case ASSIGN :
         SAY (2, "\nMatch      Do immediate assignment", Code);
         ThisStat = ASSIGN_SUBSTR;
         Header.Assign.Hdr     = StatTag;
         Header.Assign.Id      = GetWord (Pat_ptr);
         Header.Assign.Offset  = (word) (Str_ptr - Str_Beg) + Str_offset;
         goto Activate;
      case QUERY :
         SAY (2, "\nMatch      Do NOEMPTY check", Code);
         ThisStat = NOEMPTY;
         Header.Noempty.Hdr  = StatTag;
         Header.Noempty.Stub = MakeStub;
         goto Activate;
      case DO_NOT_RETURN :
         CleanMSS_and_Success;
      case INVERSE :
         SAY (2, "\nMatch      Do NOT", Code);
         ThisStat = NOT;
         Header.Not.Hdr     = StatTag;
         Header.Not.Pat_ptr = Pat_ptr;
         goto Activate;
      case LITTLE_REPEATER :
         SAY (2, "\nMatch      Do * REPEATER", Code);
         SaveStr;
         ThisStat = AS_LITTLE;
         Header.As_little.Hdr     = StatTag;
         Header.As_little.Pat_ptr = Pat_ptr;
         Header.As_little.Counter = 0;
         PushMSS;
         SkipOneItem_and_DoLittleRepeater;
         if (*Pat_ptr == DO_NOT_RETURN) Pat_ptr++;
         goto Success;
      case BIG_REPEATER :
         SAY (2, "\nMatch      Do $ REPEATER", Code);
         ThisStat = AS_MUCH;
         Header.As_much.Hdr     = StatTag;
         Header.As_much.Pat_ptr = Pat_ptr;
         Header.As_much.Counter = 0;
         SaveStr;
         goto Activate;
   }

Sleep :         /*>     Here, if it is necessary to push the statement
                        to the MSS. The statement is indicated by the
                        ThisStat variable.                              <*/
   PushMSS;

Success :       /*>     Any temporary success while matching process
                        leads here ...                                  <*/

   SetStat (ASS_ptr);
   SAY (2, "\nMatch      SUCCESS - active statement is %d ", ThisStat);
   if (ThisStat == SELECT)       /* Most frequent case                  */
      goto Process;
 /*>

   Now, the statement (ThisStat) achives its local success. Actions:

   (o)   Pop the statement header from the ASS
   (o)   Make some specific actions
   (o)   Push the header to the MSS (most frequent case)
                                                                        <*/
   if (ThisStat != AS_LITTLE) PopASS;
   switch (ThisStat)                    /* Select active statement      */
   {
      case NOT :                        /* Clean the MSS then alarm     */
         SAY (6, "\nMatch      NOT: success --> failure", Code);
         CleanMSS_and_Failure;
      case NOEMPTY :
         SAY (6, "\nMatch      NOEMPTY check ", Code);
         if (Header.Noempty.Stub == MakeStub) goto WakeUp;
         goto Sleep;
      case REPEATER :
         SAY (6, "\nMatch      FINITE REPEATER: %ld times", Header.Repeater.Counter+1);
         if ((++Header.Repeater.Counter) == Header.Repeater.Wanted) goto Sleep;
         Pat_ptr = Header.Repeater.Pat_ptr;
         goto Activate;
      case AS_MUCH :
         SAY (6, "\nMatch      $ REPEATER", Code);
         if (*Pat_ptr != DO_NOT_RETURN)
         {
            Pat_ptr = Header.As_much.Pat_ptr;
            Header.As_much.Counter++;
            SAY (6, " - %ld times OK, just one more", Header.As_much.Counter);
            goto Activate;
         }
         else     /* Disable returns if the FENCE follows the repeater  */
         {
            SAY (6, " - optimized by FENCE", Code);
            Pat_ptr = Header.As_much.Pat_ptr;
            PushASS;
            CleanMSS_and_Go;
         }
      case AS_LITTLE :
         SAY (6, "\nMatch      * REPEATER", Code);
         if (*Pat_ptr != DO_NOT_RETURN)
         {
            PopASS;
            Header.As_little.Counter++;
            SAY (6, " - %ld times OK", Header.As_little.Counter);
         }
         else
         {
            SAY (6, " - optimized by FENCE (skipped)", Code);
            Pat_ptr++;
            CleanMSS_and_ContinueLittleRepeater;
            PopASS;
         }
         SaveStr;
         goto Sleep;
      case CALL :
         SAY (6, "\nMatch      CALL: OK", Code);
         Pat_ptr = Header.Call.Pat_ptr;/* Return to parent pattern      */
         SAY (8, "\nMatch      Return to Pat_ptr = %lx (Hex)", Pat_ptr);
         goto Sleep;
      case ASSIGN_SUBSTR :
         SAY (6, "\nMatch      ASSIGN: OK", Code);
//         if (AssignVariable != NULL)
            (* AssignVariable)
            (  Header.Assign.Id,
               Header.Assign.Offset,
               (word) (Str_ptr - Str_Beg)
                  + Str_offset
/*               (word) (Str_ptr - Str_Beg)
                  + Str_offset - Header.Assign.Offset*/
            );
         ThisStat = SLEEPING_ASSIGN;
         Header.Sleeping_Assign.Hdr = StatTag;
         goto Sleep;
      case RANGE :

GeneralSuccess :

         SAY (2, "\nMatch      GENERAL SUCCESS", Code);
         *Pointer = (word) (Str_ptr - Str_Beg) + Str_offset;
         ErrorDone (MATCH_SUCCESS);
      default :
         SAY (0, "\nINTERNAL ERROR - Bad ThisStat %d [S]", ThisStat);
         goto Crash;
   }

WakeUp :        /*>     Make the statement  (ThisStat) active and
                        signal a failure for it ...                     <*/
   PushASS;

Failure :       /*>     Any temporaty failure of matching process
                        leads here ...                                  <*/

   SAY (2, "\nMatch      FAILURE ", Code);
   if (StubOnMSS)
   {/*>

      Failure  of  a  statement matching. The statement header is on the
      Active Statement Stack. The stub saves string pointer  before  the
      statement start. Actions:

      (o)   Restore string pointer
      (o)   Pop statement header from the ASS
      (o)   Switch to the statement
                                                                        <*/
      SAY (2, "- stub on MSS ", Code);
      RestoreStr;
      SAY (4, "\nMatch      Unmatch up to Str_ptr = %ld (Hex)", Str_ptr);
      SAY (4, ":'%s'", Str_ptr);
      SetStat (ASS_ptr);
      PopASS;                           /* Pop a statement from the ASS */
      switch (ThisStat)                 /* Select failed statement      */
      {
         case SELECT :                  /* Try other alternative        */
            SAY (6, "\nMatch      SELECT tries next alternative", Code);
            Pat_ptr = Header.Select.Pat_ptr;

IsAlternative :                         /* Find for next alternative    */

            if (*Pat_ptr == OR)         /* Exists ...                   */
            {
               Header.Select.Pat_ptr = ++Pat_ptr;
               goto Activate;
            }
            if (*Pat_ptr == RIGHT_BRACE) goto Failure;
            SkipOneItem_and_IsAlternative;
         case NOT :                     /* Successful matching  ...     */
            SAY (6, "\nMatch      NOT: failure --> success", Code);
            Pat_ptr = Header.Not.Pat_ptr;
            SkipOneItem_and_Success;
         case NOEMPTY :                 /* Failure matching ...         */
            SAY (6, "\nMatch      NOEMPTY failed", Code);
            goto Failure;
         case REPEATER :                /* Failure matching ...         */
            SAY (6, "\nMatch      Failure under FINITE REPEATER", Code);
            if (Header.Repeater.Counter--)
            {
               SAY (6, " - %ld times left", Header.Repeater.Counter);
               goto WakeUp;
            }
            SAY (6, " - FAILED itself", Code);
            goto Failure;
         case AS_MUCH :                 /* Always OK ...                */
            SAY (6, "\nMatch      $ REPEATER is OK", Code);
            Pat_ptr = Header.As_much.Pat_ptr;
            PushMSS;
            SkipOneItem_and_Success;
         case AS_LITTLE :               /* Can't match it               */
            SAY (6, "\nMatch      * REPEATER", Code);
            if (Header.As_little.Counter--)
            {
               SAY (6, " - %ld times left", Header.As_little.Counter);
               goto WakeUp;
            }
            else
            {
               SAY (6, " failed", Code);
               goto Failure;
            }
         case CALL :
            SAY (6, "\nMatch      CALL statement failed", Code);
	    /*>

	        The parent statement will restore its own Pat_ptr. Hence
		there is no necessity to restore Pat_ptr here.

                Pat_ptr = Header.Call.Pat_ptr;
                SAY (8, "\nMatch      Return to Pat_ptr = %lx (Hex)", Pat_ptr);

	    <*/
            goto Failure;
         case ASSIGN_SUBSTR :
            SAY (6, "\nMatch      ASSIGN failed", Code);
            goto Failure;
         default :
            SAY (0, "\nINTERNAL ERROR - Bad HDR code %o (Oct) [F]", ThisStat);
            goto Crash;
   }  }
   else
   {/*>

      Awake  a statement on the MSS stack. There are bad news for it ...
      Actions:

      (o)   Pop sleeping statement header from the MSS
      (o)   Switch to avaked statement
                                                                <*/
      SetStat (MSS_ptr);
      SAY (2, "- awake statement no.%d", ThisStat);
      PopMSS;                           /* Pop statement from the MSS   */
      switch (ThisStat)                 /* Go to the statement          */
      {
         case SELECT :
            SAY (6, "\nMatch      wake SELECT up", Code);
            goto WakeUp;
         case ELLIPSIS :
            SAY (6, "\nMatch      wake ARB (...) up", Code);
            Pat_ptr = Header.Ellipsis.Pat_ptr;
            RestoreStr;
            if (EndOfString) goto Failure;
            Str_ptr++;
            SaveStr;
            goto Sleep;
         case REPEATER :
            SAY (6, "\nMatch      wake FINITE REPEATER up", Code);
            goto WakeUp;
         case AS_MUCH :
            SAY (6, "\nMatch      wake $ REPEATER up", Code);
            SAY (6, " - %ld times left", Header.As_much.Counter - 1);
            if ((Header.As_much.Counter--) != 0) goto WakeUp;
            goto Failure;
         case AS_LITTLE :
            SAY (6, "\nMatch      wake * REPEATER up", Code);
            Pat_ptr = Header.As_little.Pat_ptr;
            RestoreStr;
            goto Activate;
         case CALL :
            SAY (6, "\nMatch      wake CALL statement up", Code);
            goto WakeUp;
         case NOEMPTY :
            SAY (6, "\nMatch      wake NOEMPTY up", Code);
            goto WakeUp;
         case FENCED_ASSIGN :
            SAY (6, "\nMatch      wake ASSIGN (fenced) up", Code);
//            if (DeAssignVariable != NULL)
               (* DeAssignVariable) (Header.Sleeping_Assign.Id);
            goto Failure;
         case SLEEPING_ASSIGN :
            SAY (6, "\nMatch      wake ASSIGN up", Code);
//            if (DeAssignVariable != NULL)
               (* DeAssignVariable) (Header.Sleeping_Assign.Id);
            ThisStat = ASSIGN_SUBSTR;
            Header.Assign.Hdr = StatTag;
            Header.Assign.Offset = (word) (Str_ptr - Str_Beg) + Str_offset;
            goto WakeUp;
         case RANGE :
            SAY (6, "\nMatch      GENERAL FAILURE (RANGE awaked)", Code);
            ErrorDone (MATCH_FAILURE);
         default :
            SAY (0, "\nINTERNAL ERROR - Bad HDR code %o (Oct) [A]", ThisStat);
            goto Crash;
   }  }

SkipOneItem :
   {/*>

      Skip one statement of the pattern. In  purpose  of  maximal  speed
      achivement this subroutine is implemented as a plain code. The Re-
      turnTo variable is used to specify return point label :

              Value     Label

                0       Success
                1       IsAlternative
                2       EndOfSelect
                3       DoLittleRepeater                        <*/

      boolean           OneMore;
      word              Nesting = 0;

      do
      {
         OneMore = FALSE;
         if ((Code = *Pat_ptr++) <= LAST_ASCII_CHARACTER)
         {
            if (Code)
            {
               while (  ((Code = *Pat_ptr++) <= LAST_ASCII_CHARACTER)
                     && Code);
               --Pat_ptr;
         }  }
         else if ((natural) Code < FIRST_STATEMENT_CODE)
         {
            if ((natural) Code <= MAX_FINITE_REPEATER)
               OneMore = TRUE;
            else if (  Code == HOLERITH
                    || Code == CASE_DEAF_HOLERITH )
            {
               Code = *Pat_ptr++;
               Pat_ptr += Code;
            }
            else if (Code == ANY_OF)
            {
               Pat_ptr += 2 + Pat_ptr [1];
            }
            else if (Code > ANY_OF)
            {
               ErrorDone (MATCH_WRONG_PATTERN_FORMAT);
         }  }
         else
         {
            switch (Code)
            {
               case LEFT_BRACE :
                  Nesting++;
                  break;
               case OR :
                  if (!OneMore && !Nesting)
                     --Pat_ptr;
                  break;
               case RIGHT_BRACE :
                  if (!OneMore && !Nesting)
                     --Pat_ptr;
                  else
                     --Nesting;
                  break;
               case EXTERNAL_PATTERN :
                  Pat_ptr += PTR_SIZE;
                  break;
               case GO_TO :
               case ASSIGN :
                  Pat_ptr += WORD_SIZE;
               case FINITE_REPEATER :
               case LITTLE_REPEATER :
               case DO_NOT_RETURN :
               case BIG_REPEATER :
               case INVERSE :
               case QUERY :
                  OneMore = TRUE;
               case ARB :
                  break;
               default :
                  ErrorDone (MATCH_WRONG_PATTERN_FORMAT);
      }  }  }
      while (OneMore || Nesting);

   }
   switch (ReturnTo)
   {
      case 0  : goto Success;
      case 1  : goto IsAlternative;
      case 2  : goto EndOfSelect;
      case 3  : goto DoLittleRepeater;
      default : SAY (0, "\nINTERNAL ERROR [R]", Code); goto Crash;
   }

CleanMSS :
   {/*>

      To avoid any return in the frame of the current alternative branch
      the MSS must be cleaned up to the stub of the last  active  state-
      ment.  In  the purpose of maximal speed achivement this subroutine
      is implemented as a plain code. The ReturnTo variable is  used  to
      specify return point label :

              Value     Label

                0       Success
                1       Failure
                2       Go
                2       ContinueLittleRepeater
                3       ContinueFence
                4       ContinueBigRepeater

              Additional bits set

              0x10      One more stub is deleted
              0x20      It is deassigned (not saved)

                                                        <*/
      register word     StubCount   = 0;
      word              AssignCount = 0;
      int               ThisStat;

      SAY (2, "\nMatch      Do FENCE", Code);
      for (;;)
      {
         if (StubOnMSS)
         {
            if (!StubCount--) break;
            KillStub;
         }
         else
         {
            SetStat (MSS_ptr);
            PopMSS;
            StubCount++;
            switch (ThisStat)
            {
               case CALL :
               case SELECT :
               case NOEMPTY :
               case ELLIPSIS :
                  continue;
               case REPEATER :
                  SAY (8, "\nMatch      Delete FINITE REPEATER", Code);
                  StubCount += Header.Repeater.Counter;
                  continue;
               case AS_MUCH :
                  SAY (8, "\nMatch      Delete $ REPEATER", Code);
                  StubCount += Header.As_much.Counter;
                  continue;
               case AS_LITTLE :
                  SAY (8, "\nMatch      Delete * REPEATER", Code);
                  StubCount += Header.As_little.Counter;
                  continue;
               case FENCED_ASSIGN :
                  --StubCount;
               case SLEEPING_ASSIGN :
                  SAY (8, "\nMatch      Deal ASSIGN", Code);
                  if (ReturnTo & 0x20)
                  {  /* Deassign                                */
                     SAY (8, "... deassigned", Code);
//                     if (DeAssignVariable != NULL)
                        (* DeAssignVariable) (Header.Sleeping_Assign.Id);
                     continue;
                  }
                  else
                  {  /* Preserve                                */
                     SAY (8, "... fenced", Code);
                     ThisStat = FENCED_ASSIGN;
                     Header.Sleeping_Assign.Hdr = StatTag;
                     PushASS;
                     AssignCount++;
                     continue;
                  }
               default :
                  SAY (0, "INTERNAL ERROR [FENCE]!", Code);
                  goto Crash;
      }  }  }
      if (ReturnTo & 0x10) KillStub;
      while (AssignCount--)     /* This is necessary for DeAssign       */
      {
         PopASS;
         PushMSS;
   }  }
   switch (ReturnTo & 0x0F)
   {
      case 0  : goto Success;
      case 1  : goto Failure;
      case 2  : goto Go;
      case 3  : goto ContinueLittleRepeater;
      default : SAY (0, "\nINTERNAL ERROR [R]", Code); goto Crash;
   }

CantReturn :

   ErrorDone (MATCH_CANNOT_RETURN);

NoMoreMemory :

   freedm ();                   /* release memory allocated by myself   */
   ErrorDone (MATCH_NO_DYNAMIC_MEMORY);

Crash :

   ErrorDone (MATCH_INTERNAL_ERROR);

}               /*>  match  <*/
/*>

   MatchKeyWord -- Compare Keyword

        String          - points to a string
        Str_End         - points to a byte followed the string
        Text            - points to a keyword. The keyword can be
                          followed by NUL, SPACE, TAB or RIGHT_ANGLE.

   Returns :

        Pointer to a character followed the keyword. A  failure  may  be
        encountered as MatchKeyWord == String.

   Matching is not case sensive. Keyword string ended by  letter,  digit
   or underline (_) must be followed by other from mentioned above char-
   acters.

<*/
#ifndef NON_ANSI

byte *  MatchKeyWord (byte * String, byte * Str_End, register byte * Text)

#else

byte *  MatchKeyWord (String, Str_End, Text)

   byte *               String;
   byte *               Str_End;
   register byte *      Text;

#endif  /*> NON_ANSI <*/
{
   register byte *      Str_ptr = String;
   register byte        Code;

   while (  (Code = (byte) tolower (*Text++)) != NUL
         && Code != RIGHT_ANGLE
         && Code != TAB
         && Code != SPACE
         )  if (  EndOfString
	       || Code != (byte) tolower (*Str_ptr++)
	       )  return (String);

   if (  !EndOfString
      && (IsUserAlpha (*(Str_ptr-1)) || isdigit (Code))
      && (IsUserAlpha (*Str_ptr)     || isdigit (Code))
      )  return String;

   return (Str_ptr);

}               /*>  MatchKeyWord  <*/

/*>

   First_pass -- First pass of the  pattern  translation.  Task  of  the
                 first pass is to build label table and to estimate size
        of pattern's internal representation. Some errors are recognized
        during this pass: brace errors, duplicate labels and so on.

        Str_ptr         - Source string pointer
        Str_End         - Points to the first byte following source string
        New_Str_ptr     - Points to the new string pointer
        Error_code      - Error code (not changed if OK)
        Error_ptr       - Points to an error location

   Returns :

        Length of pattern's internal representation

<*/
#ifndef NON_ANSI

word First_pass
(
   register byte *      Str_ptr,
   byte *               Str_End,
   byte **              New_Str_ptr,
   int  *               Error_code,
   byte **              Error_ptr
)
#else

word First_pass (Str_ptr, Str_End, New_Str_ptr, Error_code, Error_ptr)

   register byte *      Str_ptr;
   byte *               Str_End;
   byte **              New_Str_ptr;
   int  *               Error_code;
   byte **              Error_ptr;

#endif  /*> NON_ANSI <*/
{
   register byte        Code=0;
   byte *               Lex_ptr;
   byte *               Lex_End;
   word                 Pat_size     = 0;
   byte                 Zero         = 0;       /* for PushXB only      */
   boolean              AfterLiteral = FALSE;
   HeaderUnion          Header;

   MSS_mask = 0;

   while (!EndOfString)
   {
      Lex_ptr = Str_ptr;
      if (  RABBIT_EARS == (Code = *Str_ptr++)
         || APOSTROPHE  == Code
         || LEFT_ANGLE  == Code
         )
      {
         register boolean       EightBit   = FALSE;
         register byte          Quotation  = Code;
         register word          Length     = 0;
         boolean                Circumflex = FALSE;

         SAY (2, "\n1st Pass   Text constant %c", Code);
         if (Quotation == LEFT_ANGLE)
         {
            Quotation = RIGHT_ANGLE;
            EightBit  = TRUE;
         }
         while (!EndOfString && ((Code = *Str_ptr++) != Quotation))
         {
            SAY (2, "%c", Code);
            if (!Circumflex && (Code == CIRCUMFLEX))
               Circumflex = TRUE;
            else
            {
               if (Circumflex)
               {
                  //Code = DoCircumflex (Code);
                 DoCircumflex (Code);
                 Circumflex = FALSE;
               }
               Length++;
               EightBit =
                   (char)
                   (  EightBit
                   || !Code
                   || (Code > LAST_ASCII_CHARACTER)
                   );
            }
         }
         if (Code != Quotation) Error (MATCH_MISSING_QUOTATION, Str_ptr);
         SAY (2, "%c", Code);
         if (Circumflex) Length++;
         Pat_size += Length;
         if (EightBit)
         {
            AfterLiteral = FALSE;
            if (Length > MAX_HOLERITH)
               Pat_size += (  (Length + MAX_HOLERITH - 1) / MAX_HOLERITH
                           +  1
                           )
                         * 2;
            else
               Pat_size += 2;
         }
         else
         {
            if (!Length || AfterLiteral) Pat_size++;
            if (!Length)
               AfterLiteral = FALSE;
            else
               AfterLiteral = TRUE;
      }  }
      else if (Code == '{')
      {
         register int           Low  = 255;
         register int           High = 0;
         register boolean       Circumflex = FALSE;

         SAY (2, "\n1st Pass   One of the character set {", Code);
         AfterLiteral = FALSE;
         while (!EndOfString && ('}' != (Code = *Str_ptr++)))
         {
            SAY (2, "%c", Code);
            if (!Circumflex && (Code == CIRCUMFLEX))
               Circumflex = TRUE;
            else
            {
               if (Circumflex)
               {
                  //Code = DoCircumflex (Code);
                 DoCircumflex (Code);
                 Circumflex = FALSE;
               }
               if (Low  > Code) Low  = Code;
               if (High < Code) High = Code;
         }  }
         if (Code != '}') Error (MATCH_MISSING_QUOTATION, Str_ptr);
         SAY (2, "}", Code);
         if (Circumflex)
         {
            if (Low  > CIRCUMFLEX) Low  = CIRCUMFLEX;
            if (High < CIRCUMFLEX) High = CIRCUMFLEX;
         }
         Pat_size +=  3
                   +  (  Low > High
                      ?  0
                      :  1 + (High >> 3) - (Low >> 3)
                      );
      }
      else if (Code == SPACE || Code == TAB)
      {
         PassSpace;
      }
      else if (isdigit (Code))
      {
         register word  Repetition_count = 0;

         AfterLiteral = FALSE;
         SAY (2, "\n1st Pass   Finite repeater (at %lx) ", Str_ptr);
         --Str_ptr;
         do
         {
            Code = *Str_ptr++;
            if (Repetition_count <= MAX_INT_DIV_10)
               Repetition_count *= 10;
            else
               Error (MATCH_TOO_BIG_REPEATER, Lex_ptr);

            switch (Code)
            {
               case '9' : Repetition_count += 9; break;
               case '8' : Repetition_count += 8; break;
               case '7' : Repetition_count += 7; break;
               case '6' : Repetition_count += 6; break;
               case '5' : Repetition_count += 5; break;
               case '4' : Repetition_count += 4; break;
               case '3' : Repetition_count += 3; break;
               case '2' : Repetition_count += 2; break;
               case '1' : Repetition_count += 1; break;
               case '0' :                        break;
               default  : Repetition_count /= 10;
                          --Str_ptr;
                          goto CollapseRepeater;
            }
            SAY (2, "%c", Code);

         }  while (!EndOfString);

CollapseRepeater :

         if (Repetition_count > MAX_REPEATER_COUNT)
            Pat_size += 1 + WORD_SIZE;
         else
            Pat_size += 1;
         SAY (2, " (stop at %lx) ", Str_ptr);
      }
      else if (Code != LF)
      {
         SAY (2, "\n1st Pass   Keyword beginning with (%c)", Code);
         AfterLiteral = FALSE;
         switch (Code)
         {
            case '(' :
            case '[' :
               if (MSS_mask == LAST_MSS_MASK || MSS_mask == 0)
               {
                  PushStack (MSS_ptr, MSS_used, &Zero, 1);
                  MSS_mask = 1;
               }
               else
               {
                  MSS_mask <<= 1;
               }
               if (Code == '(')
               {
                  *(MSS_ptr-1) &= ~MSS_mask;
                  SAY (4,"\n           ( MSS_mask = %x (Hex)", MSS_mask);
               }
               else
               {
                  *(MSS_ptr-1) |= MSS_mask;
                  SAY (4, "\n          [ MSS_mask = %x (Hex)", MSS_mask);
               }
               Pat_size++;
               break;
            case ')' :
               PopBrace (Code);
               if (Code) Error (MATCH_BRACE_ERROR, Lex_ptr);
               Pat_size++;
               break;
            case ']' :
               PopBrace (Code);
               if (Code != 1) Error (MATCH_BRACE_ERROR, Lex_ptr);
               Pat_size += 2;
               break;
            default       :
            {
               register byte    KeyWordNo;

               for (  KeyWordNo = 0;
                      (  KeyWordNo < KEY_LIST_SIZE
                      && 0 == KeyWord ((byte *) KeyList [KeyWordNo].Name)
		      );
                      KeyWordNo++
		   );
               if (KeyWordNo < KEY_LIST_SIZE)
               {
                  if (IsUserAlpha (*Lex_ptr))
                  {
                     PassSpace;
                     if (!EndOfString)
                     switch (*Str_ptr)
                     {
                        case RIGHT_ANGLE :
                           Error (MATCH_RESERVED_KEYWORD, Lex_ptr);
                        case EQUATION    :
                           goto Assignment;
                  }  }
                  SAY (2, " recognized as %ho (Oct)",
                       KeyList [KeyWordNo].Code);
                  Pat_size++;
               }
               else
               {
                  if (IsUserAlpha (*Lex_ptr))
                  {
                     PassId;
                     if (EndOfString) goto ReferenceToLabel;
                     switch (*Str_ptr)
                     {
                        case RIGHT_ANGLE :
                        {
                           byte *       New_Str_ptr  = ++Str_ptr;
                           byte *       Old_ASS_ptr  = ASS_ptr;
                           halfword     Old_ASS_used = ASS_used;

                           if (  GetExternalPattern != NULL
                              && (* GetExternalPattern)
                                 (  (char *) Lex_ptr,
                                    Lex_End - Lex_ptr
                                 ) != NULL
                              )  Error (MATCH_RESERVED_KEYWORD, Lex_ptr);

                           SAY (2, " defines a label", Code);
                           while (!EmptyASS)
                           {
                              PopLabel;
                              if (KeyWord (Header.Label.Name))
                                 Error (MATCH_DUPLICATE_LABEL, Lex_ptr);
                           }
                           ASS_ptr  = Old_ASS_ptr;
                           ASS_used = Old_ASS_used;
                           Header.Label.Name = Lex_ptr;
                           Header.Label.Destination_offset = Pat_size;
                           PushLabel;
                           Str_ptr = New_Str_ptr;
                        }  break;
                        case EQUATION :
Assignment :               Str_ptr++;
                        default :
ReferenceToLabel :         Pat_size += 1 + WORD_SIZE;
                           SAY (2, " is label or equation", Code);
                  }  }
                  else
                  {
                     PopBrace (Code);
                     if (Code != 2)
                        Error (MATCH_UNRECOGNIZED_CHARACTER, Lex_ptr);

                     *New_Str_ptr = Str_ptr;
                     return (++Pat_size);
   }  }  }  }  }  }

   PopBrace (Code);
   if (Code != 2) Error (MATCH_MISSING_RIGHT_BRACE, Str_ptr);
   *New_Str_ptr = Str_ptr;
   return (++Pat_size);

NoMoreMemory :

   freedm ();                   /* release memory allocated by myself   */
   Error (MATCH_NO_DYNAMIC_MEMORY, Str_ptr);

}               /*>  First_pass  <*/

/*>

   CloseStatement -- Pop the MSS until a select statement

        ResultStatus    - EMPTY_ATOM if a subpattern matches null string,
                          NON_EMPTY_ATOM - otherwise

   Returns :

      EMPTY_ATOM        if there is possibility to match null string
      NON_EMPTY_ATOM    if there is no possibility to match null
      ERROR_STATUS      if $ or * repeat a pattern which can match null

<*/
#ifndef NON_ANSI

byte    CloseStatement (register byte ResultStatus)

#else

byte    CloseStatement (ResultStatus)

   register byte        ResultStatus;

#endif  /*> NON_ANSI <*/
{
   byte         EmptyStatus;

   while (GetHdr (MSS_ptr) & NON_SELECT_MASK)
   {
      PopStatement;
      switch (EmptyStatus)
      {
         case INDEFINITE :
            if (ResultStatus == EMPTY_ATOM) return (ERROR_STATUS);
            ResultStatus = EMPTY_ATOM;
         case SON_WILL_SAY :
            break;
         default :
            ResultStatus = (byte) (EmptyStatus >> NON_SELECT_MASK);
   }  }
   *(MSS_ptr-1) |= ResultStatus;
   return (ResultStatus);

}               /*>  CloseStatement  <*/
/*>

   Second_pass -- Second pass of the pattern translation.

        Str_ptr         - Source string pointer
        Str_End         - Points to the first byte following source string
        Pat_ptr         - Points to a buffer for the pattern allocation
        Error_code      - Error code (not changed if OK)
        Error_ptr       - Points to an error location

<*/
#ifndef NON_ANSI

int Second_pass
(
   register byte *      Str_ptr,
   byte          *      Str_End,
   register byte *      Pat_ptr,
   int           *      Error_code,
   byte          **     Error_ptr
)
#else

int Second_pass (Str_ptr, Str_End, Pat_ptr, Error_code, Error_ptr)

   register byte *      Str_ptr;
   byte          *      Str_End;
   register byte *      Pat_ptr;
   int           *      Error_code;
   byte          **     Error_ptr;

#endif  /*> NON_ANSI <*/

{
   register byte        Code;
   byte          *      Lex_ptr;
   byte          *      Lex_End;
   byte          *      Pat_Beg = Pat_ptr;
   byte                 EmptyStatus = SELECT_BRANCH;
   boolean              AfterLiteral = FALSE;
   HeaderUnion          Header;

   PushStatement;
   while (!EndOfString)
   {
      Lex_ptr = Str_ptr;
      if ((Code = *Str_ptr++) != SPACE && Code != TAB && Code != LF)
      {
         EmptyStatus = NON_EMPTY_ATOM;
         if (  Code == RABBIT_EARS
            || Code == APOSTROPHE
            || Code == LEFT_ANGLE
            )
         {
            register boolean    EightBit   = FALSE;
            register byte       Quotation  = Code;
            register word       Length     = 0;
            word                Counter;
            boolean             Broken;
            boolean             Circumflex = FALSE;

            SAY (4, "\n2nd Pass   Text constant %c", Code);
            if (Quotation == LEFT_ANGLE)
            {
               Quotation = RIGHT_ANGLE;
               EightBit  = TRUE;
            }
            while (!EndOfString && ((Code = *Str_ptr++) != Quotation))
            {
               SAY (4, "%c", Code);
               if (!Circumflex && (Code == CIRCUMFLEX))
                  Circumflex = TRUE;
               else
               {
                  if (Circumflex)
                  {
                     //Code = DoCircumflex (Code);
                    DoCircumflex (Code);
                    Circumflex = FALSE;
                  }
                  if (Quotation == RIGHT_ANGLE) Code = (byte) tolower (Code);
                  Length++;
                  EightBit =
		     (byte)
		     (  EightBit
		     || !Code
		     || (Code > LAST_ASCII_CHARACTER)
		     );
            }  }
            SAY (4, "%c", Code);
            if (Circumflex) Length++;
            if (!Length) EmptyStatus = EMPTY_ATOM;
            Str_ptr = ++Lex_ptr;
            Broken = (boolean) (Length > MAX_HOLERITH);
            SAY (2, "\n2nd Pass   << ", Code);
            if (Broken)
            {
               *Pat_ptr++ = LEFT_BRACE;
               SAY (2, "(", Code);
            }
            if (Broken || EightBit)
            {
               AfterLiteral = FALSE;
            }
            else
            {
               if (!Length || AfterLiteral)
               {
                  *Pat_ptr++ = NOOP;
                  SAY (6," [NOOP to separate literals]", Code);
               }
               if (!Length)
                  AfterLiteral = FALSE;
               else
                  AfterLiteral = TRUE;
            }
            Counter = Length;
            while (Counter)
            {
               Length = Counter;
               if (Broken || EightBit)
               {
                  if (Length > MAX_HOLERITH) Length = MAX_HOLERITH;
                  if (Quotation == RIGHT_ANGLE)
                  {
                     *Pat_ptr++ = CASE_DEAF_HOLERITH;
                     SAY (2, " CASE DEAF [%d] ", Length - 1);
                  }
                  else
                  {
                     *Pat_ptr++ = HOLERITH;
                     SAY (2, " HOLERITH [%d] ", Length - 1);
                  }
                  *Pat_ptr++ = (byte) (Length - 1);
               }
               Counter -= Length;
               Circumflex = FALSE;
               while (Length)
               {
                  Code = *Str_ptr++;
                  if (!Circumflex && (Code == CIRCUMFLEX))
                     Circumflex = TRUE;
                  else
                  {
                     --Length;
                     if (Circumflex)
                     {
                        *Pat_ptr++ = DoCircumflex (Code);
                        SAY (2, "%c", DoCircumflex (Code));
                        Circumflex = FALSE;
                     }
                     else
                     {
                        if (Quotation == RIGHT_ANGLE)
			{
			   Code = (byte) tolower (Code);
			}
                        *Pat_ptr++ = Code;
                        SAY (2, "%c", Code);
            }  }  }  }
            if (EmptyStatus != EMPTY_ATOM && *(Str_ptr-1) == Quotation)
            {
               *(Pat_ptr-1) = CIRCUMFLEX; /* End cicumflex means itself */
               SAY (2, "%c", CIRCUMFLEX);
            }
            else
            {
               Str_ptr++;               /* Pass closing quotation marks */
            }
            if (Broken)
            {
               *Pat_ptr++ = RIGHT_BRACE;
               SAY (2, ")", Code);
         }  }
         else if (Code == '{')
         {
            register int        Low  = 255;
            register int        High = 32;
            register boolean    Circumflex = FALSE;
            static byte         Map [32];

            SAY (2, "\n2nd Pass   One of the character set {", Code);
            while (0 != --High) Map [High] = 0;
            AfterLiteral = FALSE;
            EmptyStatus  = RETURN_NON_EMPTY;
            while (!EndOfString && ('}' != (Code = *Str_ptr++)))
            {
               SAY (2, "%c", Code);
               if (!Circumflex && (Code == CIRCUMFLEX))
                  Circumflex = TRUE;
               else
               {
                  if (Circumflex)
                  {
                     //Code = DoCircumflex (Code);
                    DoCircumflex (Code);
                    Circumflex = FALSE;
                  }
                  if (Low  > Code) Low  = Code;
                  if (High < Code) High = Code;
                  Map [Code >> 3] |= PowerOf2 [Code & 0x07];
            }  }
            SAY (2, "} =", Code);
            if (Circumflex)
            {
               if (Low  > CIRCUMFLEX) Low  = CIRCUMFLEX;
               if (High < CIRCUMFLEX) High = CIRCUMFLEX;
               Map [CIRCUMFLEX >> 3] |= PowerOf2 [CIRCUMFLEX & 0x07];
            }
            *Pat_ptr++ = ANY_OF;
            if (Low > High)
            {
               *Pat_ptr++ = 0;
               *Pat_ptr++ = 0;
               SAY (2, " nothing", Code);
            }
            else
            {
               *Pat_ptr++ = (byte) (Low >>= 3);
               *Pat_ptr++ = (byte) (High = 1 + (High >> 3) - Low);
               while (High--)
               {
                  SAY (2, " %x", Map [Low]);
                  *Pat_ptr++ = Map [Low++];
         }  }  }
         else if (isdigit (Code))
         {
            register word       Repetition_count = 0;

            AfterLiteral = FALSE;
            --Str_ptr;
            do
            {
               Code = *Str_ptr++;
               Repetition_count *= 10;
               switch (Code)
               {
                  case '9' : Repetition_count += 9; break;
                  case '8' : Repetition_count += 8; break;
                  case '7' : Repetition_count += 7; break;
                  case '6' : Repetition_count += 6; break;
                  case '5' : Repetition_count += 5; break;
                  case '4' : Repetition_count += 4; break;
                  case '3' : Repetition_count += 3; break;
                  case '2' : Repetition_count += 2; break;
                  case '1' : Repetition_count += 1; break;
                  case '0' :                        break;
                  default  : Repetition_count /= 10;
                             --Str_ptr;
                             goto CollapseRepeater;
               }
            }  while (!EndOfString);

CollapseRepeater :

            switch (Repetition_count)
            {
               case  0 : EmptyStatus = RETURN_EMPTY;    break;
               default : EmptyStatus = SON_WILL_SAY;
            }
            if (Repetition_count <= MAX_REPEATER_COUNT)
            {
               *Pat_ptr++ = (byte) (Repetition_count + LAST_ASCII_CHARACTER + 1);
               SAY (2, "\n2nd Pass   << repeat %d (small)", Repetition_count);
            }
            else
            {
               *Pat_ptr++ = FINITE_REPEATER;
               Take.Int = Repetition_count;
               PutInt (Pat_ptr);
               SAY (2, "\n2nd Pass   << repeat %d (big)", Repetition_count);
         }  }
         else
         {
            SAY (2, "\n2nd Pass   Keyword beginning with (%c)", Code);
            AfterLiteral = FALSE;
            switch (Code)
            {
               case ']' : CloseBranch(EMPTY_ATOM); *Pat_ptr++ = OR;
               case ')' : CloseBranch(EMPTY_ATOM); *Pat_ptr++ = RIGHT_BRACE;
                  PopStatement;
                  EmptyStatus =
		     (byte)
		     ((EmptyStatus & ~SELECT_MASK) >> SELECT_SHIFT);
                  break;
               case '(' :
               case '[' : *Pat_ptr++ = LEFT_BRACE;
                  EmptyStatus = SELECT_BRANCH;
                  break;
               default  :
               {
                  register byte KeyWordNo;

                  for (  KeyWordNo = 0;
                         (  KeyWordNo < KEY_LIST_SIZE
                         && 0 == KeyWord ((byte *) KeyList [KeyWordNo].Name)
			 );
                         KeyWordNo++
		      );
                  if (KeyWordNo < KEY_LIST_SIZE)
                  {
                     if (IsUserAlpha (*Lex_ptr))
                     {
                        PassSpace;
                        if (!EndOfString && *Str_ptr == EQUATION)
                           goto Assignment;
                     }
                     SAY (2, " recognized as %ho (Oct)", KeyList [KeyWordNo].Code);
                     *Pat_ptr++ = KeyList [KeyWordNo].Code;
                     switch (KeyList [KeyWordNo].Code)
                     {
                        case OR :
                           CloseBranch (EMPTY_ATOM);
                           continue;
                        case ELLIPSIS :
                        case END_OF_STRING :
                        case END_OF_KEYWORD :
                           EmptyStatus = EMPTY_ATOM;
                           break;
                        case BIG_REPEATER :
                        case LITTLE_REPEATER :
                           EmptyStatus = INDEFINITE;
                           break;
                        case QUERY :
                           EmptyStatus = RETURN_NON_EMPTY;
                           break;
                        case INVERSE :
                           EmptyStatus = RETURN_EMPTY;
                           break;
                  }  }
                  else
                  {
                     PassId;
                     if (EndOfString) goto ReferenceToLabel;
                     switch (*Str_ptr)
                     {
                        case RIGHT_ANGLE :
                           Str_ptr++;
                           SAY (2, " is label declaration", Code);
                           continue;
                        case EQUATION :
Assignment :            {
                           word         Id;

                           Str_ptr++;
                           SAY (2, " is immediate assignment", Code);
                           if ( ( Id = (GetVariableId)
                                        (  (char *) Lex_ptr,
                                           Lex_End - Lex_ptr
                                        ) == 0 )
                              )
                           {
                              Error (MATCH_UNDEFINED_VARIABLE, Lex_ptr);
                           }
                           *Pat_ptr++ = ASSIGN;
                           PutWord (Pat_ptr, Id);
                           EmptyStatus = SON_WILL_SAY;
                        }  break;
                        default :
ReferenceToLabel :      {
                           byte *         New_Str_ptr  = Str_ptr;
                           byte *         Old_ASS_ptr  = ASS_ptr;
                           halfword       Old_ASS_used = ASS_used;

                           SAY (2, " is reference to a label", Code);
Find_for_label :
                           if (EmptyASS)
                           {
                              char *   Pattern;

                              SAY (2, " EXTERNAL one!", Code);
                              if (  GetExternalPattern == NULL
                                 || (Pattern = (* GetExternalPattern)
                                               (  (char *) Lex_ptr,
                                                  Lex_End - Lex_ptr)
                                               ) == NULL
                                 )
                              {
                                 Error (MATCH_UNRECOGNIZED_KEYWORD, Lex_ptr);
                              }
                              *Pat_ptr++ = EXTERNAL_PATTERN;
                              PutPtr (Pat_ptr, (byte *) Pattern);
                           }
                           else
                           {
                              PopLabel;
                              if (KeyWord (Header.Label.Name))
                              {
                                  *Pat_ptr++ = GO_TO;
                                  Take.Int =
                                     (  (long) (Pat_Beg - Pat_ptr)
                                     -  WORD_SIZE
                                     +  Header.Label.Destination_offset
                                     );
                                  SAY (4, " (rel.offs. = %d )", Take.Int);
                                  PutInt (Pat_ptr);
                              }
                              else
                                 goto Find_for_label;
                           }
                           ASS_ptr  = Old_ASS_ptr;
                           ASS_used = Old_ASS_used;
                           Str_ptr  = New_Str_ptr;
         }  }  }  }  }  }
         switch (EmptyStatus)
         {
            case EMPTY_ATOM :
            case NON_EMPTY_ATOM :
               FindForSelect(EmptyStatus);
               break;
            default :
               PushStatement;
   }  }  }
   *Pat_ptr = RIGHT_BRACE;
   return (0);

NoMoreMemory :

   freedm ();                   /* release memory allocated by myself   */
   Error (MATCH_NO_DYNAMIC_MEMORY, Str_ptr);

IndefiniteLoop :

   Error (MATCH_POSSIBLE_INDEFINITE_LOOP, --Str_ptr);

}               /*>  Second_pass  <*/
/*>

   patran -- pattern translator

        Length  - of a character string containing a pattern to be translated
        String  - points to the string (maybe not NUL terminated)
        Pointer - number of starting character (0..Length)
        Pattern - address of a buffer for translated pattern
        Size    - its size

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

int patran
(
   int          Length,
   const char * String,
   int        * Pointer,
   char       * Pattern,
   int        * Size
)
#else

int patran (Length, String, Pointer, Pattern, Size)

   int          Length;
   char       * String;
   int        * Pointer;
   char       * Pattern;
   int        * Size;

#endif  /*> NON_ANSI <*/
{
   register byte *      Str_ptr = (byte *) &String [*Pointer];
   register byte *      Pat_ptr = (byte *) Pattern;
   byte          *      Str_End = (byte *) &String [Length];
   byte          *      Str_Beg = (byte *) String;
   word                 Pat_size;
   byte          *      Error_ptr;
   halfword             Old_MSS_mask = MSS_mask;
   StackData;

   CheckString;
   SaveStack;
   MatchError = MATCH_SUCCESS;
   Pat_size = First_pass (Str_ptr, Str_End, &Str_End, &MatchError, &Error_ptr);
   MSS_mask = Old_MSS_mask;
   if (MatchError != MATCH_SUCCESS) goto FatalError;
   if (Pat_size > MAX_STRING_LENGTH) ErrorDone (MATCH_TOO_LARGE_STRING);
   *Size = (int) Pat_size;
   *Pointer = Str_End - Str_Beg;
   if ((int) Pat_size > *Size)
   {
      *Size = 0;
      ErrorDone (MATCH_TOO_LARGE_PATTERN);
   }
   MSS_ptr  = Old_MSS_ptr;              /* Clean MSS after 1st Pass     */
   MSS_used = Old_MSS_used;
   (void) Second_pass (Str_ptr, Str_End, Pat_ptr, &MatchError, &Error_ptr);
   if (MatchError != MATCH_SUCCESS) goto FatalError;
   Done (MATCH_SUCCESS);

FatalError :

   *Size = 0;
   *Pointer = Error_ptr - Str_Beg;
   Done (MatchError);

}               /*>  patran  <*/
/*>

   patmaker -- pattern constructor

        String  - points to the NUL terminated string

   Returns :

        Pointer to translated pattern or NULL

   The String is translated and a memory block is requested to  allocate
   the  pattern  body.  Any syntax or other error leads to the NULL as a
   result. Note that whole String must be recognized as a pattern. Error
   code  is returned to the MatchError variable.  You should member that
   error may be resulted from malloc's failure.

<*/
#ifndef NON_ANSI

extern const char * patmaker (const char * String)

#else

extern char * patmaker (String)

   char *       String;

#endif  /*> NON_ANSI <*/
{
   byte *       Pat_ptr;
   byte *       Str_ptr      = (byte *) String;
   byte *       Str_End      = (byte *) String;
   word         Pat_size;
   byte *       Error_ptr;
   halfword     Old_MSS_mask = MSS_mask;
   StackData;

   SaveStack;
   MatchError = MATCH_SUCCESS;
   Str_End += strlen (String);
   Pat_size = First_pass (Str_ptr, Str_End, &Str_ptr, &MatchError, &Error_ptr);
   MSS_mask = Old_MSS_mask;
   if (MatchError != MATCH_SUCCESS) goto Stop;
   if (Str_ptr != Str_End)
   {
      MatchError = MATCH_UNRECOGNIZED_KEYWORD;
      goto Stop;
   }
   if (Pat_size > MAX_STRING_LENGTH)
   {
      MatchError = MATCH_TOO_LARGE_PATTERN;
      goto Stop;
   }
   if (NULL == (Pat_ptr = (byte *) GC_malloc (Pat_size)))
   {
      MatchError = MATCH_NO_DYNAMIC_MEMORY;
      goto Stop;
   }
   MSS_ptr  = Old_MSS_ptr;              /* Clean MSS after 1st Pass     */
   MSS_used = Old_MSS_used;
   (void) Second_pass ((byte *) String, Str_End, Pat_ptr, &MatchError, &Error_ptr);
   if (MatchError == MATCH_SUCCESS) Done ((char *) Pat_ptr);
   GC_free ((char *) Pat_ptr);

Stop :

   Done (NULL);

}               /*>  patmaker  <*/
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
const char * patinit
(
   int		Length,
   const char *	String,
   int        *	Pointer
)
#else
char * patinit (Length, String, Pointer)
   int		Length;
   char       *	String;
   int        *	Pointer;
#endif	/*> NON_ANSI <*/
{
   byte      *	Str_ptr = (byte *) &String [*Pointer];
   byte      *	Pat_ptr;
   byte      *	Str_End = (byte *) &String [Length];
   byte      *	Str_Beg = (byte *) String;
   word         Pat_size;
   byte      *	Error_ptr;
   halfword	Old_MSS_mask = MSS_mask;
   StackData;

   CheckString;
   SaveStack;
   MatchError = MATCH_SUCCESS;
   Pat_size =
      First_pass
      (  Str_ptr,
         Str_End,
         &Str_End,
         &MatchError,
         &Error_ptr
      );
   MSS_mask = Old_MSS_mask;
   if (MatchError != MATCH_SUCCESS) goto FatalError;
   if (  Pat_size > MAX_STRING_LENGTH
      || 0 == AllocatePattern
      || (  0
         == (  Pat_ptr
             = (byte *)
	       (*AllocatePattern) ((unsigned long) Pat_size)
      )  )  )
   {  /* No memory for the pattern				*/
      MatchError = MATCH_TOO_LARGE_PATTERN;
      Error_ptr = Str_End;
      goto FatalError;
   }
   MSS_ptr  = Old_MSS_ptr;	/* Clean MSS after 1st Pass	*/
   MSS_used = Old_MSS_used;
   *Pointer = Str_End - Str_Beg;
   (void) Second_pass
   (  Str_ptr,
      Str_End,
      Pat_ptr,
      &MatchError,
      &Error_ptr
   );
   if (MatchError != MATCH_SUCCESS) goto FatalError;
   Done ((char *) Pat_ptr);

FatalError :

   *Pointer = Error_ptr - Str_Beg;
   Done (0);

}               /*>  patinit  <*/
