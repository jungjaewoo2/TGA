* sparse.f:  Fortran sparse matrix subroutines
*
* Copyright (C) 1994-95  K. Scott Hunziker.
* Copyright (C) 1990-94  The Boeing Company.
*
* See the file COPYING (below) for license, warranty, and permission details.
*
*
* Alki is free software.  You can redistribute it and/or modify it
* under the terms of the GNU General Public License as published by
* the Free Software Foundation; either version 2 of the License, or
* (at your option) any later version.
*
* Alki is distributed in the hope that it will be useful, but WITHOUT
* ANY WARRANTY; without even the implied warranty of MERCHANTABILITY
* or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public
* License for more details.
*
* You should have received a copy of the GNU General Public License
* along with Alki; see the file LICENSE.  If not, write to the Free
* Software Foundation, Inc., 675 Mass Ave, Cambridge, MA 02139, USA.
*
* The copyright to major portions of Alki belongs to The Boeing
* Company.  The following permission notice and warranty disclaimer
* pertain to those portions of the code:
*
*    Permission to use, copy, modify, and distribute this software
*    and its documentation for any purpose is hereby granted,
*    provided that the above copyright notice appears in all copies,
*    that both the copyright notice and this permission notice and
*    warranty disclaimer appear in supporting documentation, and that
*    the names of Boeing or any of its entities not be used in
*    advertising or publicity pertaining to distribution of the
*    software without specific, written, prior permission.
*
*    BOEING DISCLAIMS ALL WARRANTIES WITH REGARD TO THIS SOFTWARE,
*    INCLUDING ALL IMPLIED WARRANTIES OF MERCHANTABILITY, FITNESS,
*    AND NONINFRINGEMENT.  IN NO EVENT SHALL BOEING BE LIABLE FOR ANY
*    SPECIAL, DIRECT, INDIRECT, OR CONSEQUENTIAL DAMAGES OR ANY
*    DAMAGES WHATSOEVER RESULTING FROM LOSS OF USE, DATA, OR PROFITS,
*    WHETHER IN AN ACTION OF CONTRACT, NEGLIGENCE, OR OTHER TORTIOUS
*    ACTION, ARISING OUT OF OR IN CONNECTION WITH THE USE OR
*    PERFORMANCE OF THIS SOFTWARE.
*
***********************************************************************

      SUBROUTINE DGSADD( AI, AJ, AN, BI, BJ, BN, N, M,
     $     CI, CJ, CN, W )
      INTEGER AI(*), AJ(*), BI(*), BJ(*), N, M
      INTEGER CI(*), CJ(*)
      DOUBLE PRECISION AN(*), BN(*), CN(*), W(M)

*     DGSADD adds two double precision matrices in the full, unordered,
*     sparse form.  The CI and CJ vectors should be previously set by a
*     call to XGSADD.  This routine does not change the order of CI and
*     CJ, so one can call XGSORD to order them before calling this
*     routine.

*     Input:   AI, AJ, AN       The first matrix.
*              BI, BJ, BN       The second matrix.
*              CI, CJ           The structure of the result.
*              N                The number of rows.
*              M                The number of columns.

*     Output:  CN               The resulting matrix.
*              W                Workspace of length M.

      INTEGER I, J, CI1, CI2

*     Loop over all the rows.

      DO 100 I = 1, N

         CI1 = CI(I)
         CI2 = CI(I+1) - 1

*        Any nonzeros in this row?

         IF ( CI1 .GT. CI2 ) GO TO 100

*        Set result columns to zero.

         DO 10 J = CI1, CI2
            W(CJ(J)) = 0.0D0
 10      CONTINUE

*        Assign nonzeros from left matrix.

         DO 20 J = AI(I), AI(I+1) - 1
            W(AJ(J)) = AN(J)
 20      CONTINUE

*        Add in nonzeros from right matrix.

         DO 30 J = BI(I), BI(I+1) - 1
            W(BJ(J)) = W(BJ(J)) + BN(J)
 30      CONTINUE

*        Collect the results.

         DO 40 J = CI1, CI2
            CN(J) = W(CJ(J))
 40      CONTINUE

 100  CONTINUE

      RETURN
      END

***********************************************************************

      SUBROUTINE DGSMUL( AI, AJ, AN, BI, BJ, BN, NA, MB,
     $     CI, CJ, CN, W )
      INTEGER AI(*), AJ(*), BI(*), BJ(*), NA, MB
      INTEGER CI(*), CJ(*)
      DOUBLE PRECISION AN(*), BN(*), CN(*), W(MB)

*     DGSMUL multiplies two matrices in the full, unordered, sparse
*     form.  The structure of the result must have already been formed,
*     likely by calling XGSMUL.  This routine does not change the order
*     of CI and CJ, so one can call XGSORD to order them before calling
*     this routine.

*     Input:   AI, AJ, AN       The first matrix.
*              BI, BJ, BN       The second matrix.
*              CI, CJ           The results matrix structure.
*              NA               The number of rows of A.
*              MB               The number of columns of B.

*     Output:  CN               The resulting matrix.
*              W                Workspace of length MB.

      INTEGER I, J, K, CI1, CI2

*     Loop over the rows of A.

      DO 100 I = 1, NA

         CI1 = CI(I)
         CI2 = CI(I+1) - 1

*        Any nonzeros in this row?

         IF ( CI1 .GT. CI2 ) GO TO 100

*        Set result columns to zero.

         DO 10 J = CI1, CI2
            W(CJ(J)) = 0.0D0
 10      CONTINUE

*        Step through the nonzeros in this row of A.

         DO 30 J = AI(I), AI(I+1) - 1

*           Work through the nonzeros in the corresponding row of B.

            DO 20 K = BI(AJ(J)), BI(AJ(J)+1) - 1
               W(BJ(K)) = W(BJ(K)) + AN(J)*BN(K)
 20         CONTINUE

 30      CONTINUE

*        Collect the results.

         DO 40 J = CI1, CI2
            CN(J) = W(CJ(J))
 40      CONTINUE

 100  CONTINUE

      RETURN
      END

***********************************************************************

      SUBROUTINE DGSSUB( AI, AJ, AN, BI, BJ, BN, N, M,
     $     CI, CJ, CN, W )
      INTEGER AI(*), AJ(*), BI(*), BJ(*), N, M
      INTEGER CI(*), CJ(*)
      DOUBLE PRECISION AN(*), BN(*), CN(*), W(M)

*     DGSSUB subtracts two double precision matrices in the full,
*     unordered, sparse form.  The CI and CJ vectors should be
*     previously set by a call to XGSADD.  This routine does not change
*     the order of CI and CJ, so one can call XGSORD to order them
*     before calling this routine.

*     Input:   AI, AJ, AN       The first matrix.
*              BI, BJ, BN       The second matrix.
*              CI, CJ           The structure of the result.
*              N                The number of rows.
*              M                The number of columns.

*     Output:  CN               The resulting matrix.
*              W                Workspace of length M.

      INTEGER I, J, CI1, CI2

*     Loop over all the rows.

      DO 100 I = 1, N

         CI1 = CI(I)
         CI2 = CI(I+1) - 1

*        Any nonzeros in this row?

         IF ( CI1 .GT. CI2 ) GO TO 100

*        Set result columns to zero.

         DO 10 J = CI1, CI2
            W(CJ(J)) = 0.0D0
 10      CONTINUE

*        Assign nonzeros from left matrix.

         DO 20 J = AI(I), AI(I+1) - 1
            W(AJ(J)) = AN(J)
 20      CONTINUE

*        Subtract nonzeros from right matrix.

         DO 30 J = BI(I), BI(I+1) - 1
            W(BJ(J)) = W(BJ(J)) - BN(J)
 30      CONTINUE

*        Collect the results.

         DO 40 J = CI1, CI2
            CN(J) = W(CJ(J))
 40      CONTINUE

 100  CONTINUE

      RETURN
      END

***********************************************************************

      SUBROUTINE DGSTRN( AI, AJ, AN, N, M, ATI, ATJ, ATN )
      INTEGER AI(*), AJ(*), ATI(*), ATJ(*), N, M
      DOUBLE PRECISION AN(*), ATN(*)

*     DGSTRN transposes a double precision matrix in the full,
*     unordered, sparse form.

*     Input:   AI, AJ, AN       The original matrix.
*              N                The number of rows.
*              M                The number of columns.

*     Output:  ATI, ATJ, ATN    The transposed matrix.

      INTEGER I, J, K, JJ

      DO 10 I = 2, M+1
         ATI(I) = 0
 10   CONTINUE

      DO 20 I = 1, AI(N+1) - 1
         J = AJ(I) + 2
         IF ( J .LE. M+1 ) ATI(J) = ATI(J) + 1
 20   CONTINUE

      ATI(1) = 1
      ATI(2) = 1

      DO 30 I = 3, M+1
         ATI(I) = ATI(I) + ATI(I-1)
 30   CONTINUE

      DO 50 I = 1, N

         DO 40 J = AI(I), AI(I+1) - 1

            JJ = AJ(J) + 1
            K = ATI(JJ)

            ATJ(K) = I
            ATN(K) = AN(J)
            ATI(JJ) = K + 1

 40      CONTINUE

 50   CONTINUE

      RETURN
      END

***********************************************************************

      SUBROUTINE XGSADD( AI, AJ, BI, BJ, N, M, CI, CJ, IW )
      INTEGER AI(*), AJ(*), BI(*), BJ(*), N, M
      INTEGER CI(*), CJ(*), IW(M)

* XGSADD symbolicly adds two matrices in the full, unordered sparse
* form.

*     Input:   AI, AJ  The first matrix structure.
*              BI, BJ  The second matrix structure.
*              N       The number of rows.
*              M       The number of columns.

*     Output:  CI, CJ  The resulting matrix structure.
*              IW      Integer workspace of length M.

      INTEGER I, J, K

      K = 1

      DO 10 I = 1, M
         IW(I) = 0
 10   CONTINUE

      DO 40 I = 1, N

         CI(I) = K

         DO 20 J = AI(I), AI(I+1) - 1
            CJ(K) = AJ(J)
            IW(AJ(J)) = I
            K = K + 1
 20      CONTINUE

         DO 30 J = BI(I), BI(I+1) - 1
            IF ( IW(BJ(J)) .EQ. I ) GO TO 30
            CJ(K) = BJ(J)
            K = K + 1
 30      CONTINUE

 40   CONTINUE

      CI(N+1) = K

      RETURN
      END

***********************************************************************

      SUBROUTINE XGSMUL( AI, AJ, BI, BJ, NA, MB, CI, CJ, MAXCJ, IW )
      INTEGER AI(*), AJ(*), BI(*), BJ(*), NA, MB
      INTEGER CI(*), CJ(MAXCJ), IW(MB)

*     XGSMUL symbolicly multiplies two matrices in full, unordered
*     sparse form.  If CJ is not dimensioned large enough, we return
*     with CI(1)=0.

*     Input:   AI, AJ  The first matrix structure.
*              BI, BJ  The second matrix structure.
*              NA      The number of rows of A.
*              MB      The number of columns of B.

*     Output:  CI, CJ  The resulting matrix structure.
*              MAXCJ   The dimensioned length of CJ.
*              IW      Integer workspace of length MB.

      INTEGER I, J, K, JJ
*
      JJ = 1

*     Initialize

      DO 10 I = 1, MB
         IW(I) = 0
 10   CONTINUE

*     Loop over rows of left matrix.

      DO 40 I = 1, NA

         CI(I) = JJ

         DO 30 J = AI(I), AI(I+1) - 1

            DO 20 K = BI(AJ(J)), BI(AJ(J)+1) - 1

               IF ( IW(BJ(K)) .EQ. I ) GO TO 20

               IF ( JJ .GT. MAXCJ ) THEN
                  CI(1) = 0
                  GO TO 50
               ENDIF

               CJ(JJ) = BJ(K)
               JJ = JJ + 1
               IW(BJ(K)) = I

 20         CONTINUE
 30      CONTINUE
 40   CONTINUE

      CI(NA+1) = JJ

 50   CONTINUE
      RETURN
      END

***********************************************************************

      SUBROUTINE XGSTRN( AI, AJ, N, M, ATI, ATJ )
      INTEGER AI(*), AJ(*), ATI(*), ATJ(*), N, M

*     XGSTRN symbolically transposes a matrix in the full, unordered,
*     sparse form.  In the process, the rows of the matrix are ordered,
*     so that two calls to this routine will order the matrix.

*     Input:   AI, AJ           The original matrix structure.
*              N                The number of rows.
*              M                The number of columns.

*     Output:  ATI, ATJ         The transposed matrix structure.

      INTEGER I, J

      DO 10 I = 2, M+1
         ATI(I) = 0
 10   CONTINUE

      DO 20 I = 1, AI(N+1) - 1
         J = AJ(I) + 2
         IF ( J .LE. M+1 ) ATI(J) = ATI(J) + 1
 20   CONTINUE

      ATI(1) = 1
      ATI(2) = 1

      DO 30 I = 3, M+1
         ATI(I) = ATI(I) + ATI(I-1)
 30   CONTINUE

      DO 50 I = 1, N

         DO 40 J = AI(I), AI(I+1) - 1
            ATJ( ATI( AJ(J)+1 ) ) = I
            ATI( AJ(J)+1 ) = ATI( AJ(J)+1 ) + 1
 40      CONTINUE

 50   CONTINUE

      RETURN
      END

***********************************************************************

      SUBROUTINE ZGSADD( AI, AJ, AN, BI, BJ, BN, N, M,
     $     CI, CJ, CN, W )
      INTEGER AI(*), AJ(*), BI(*), BJ(*), N, M
      INTEGER CI(*), CJ(*)
      DOUBLE PRECISION AN(*), BN(*), CN(*), W(*)

*     ZGSADD adds matrices in the full, unordered, sparse form.  The CI
*     and CJ vectors should be previously set by a call to XGSADD.  This
*     routine does not change the order of CI and CJ, so one could call
*     XGSORD to order them before calling this routine.

*     Double precision complex data is stored as two consecutive double
*     precision numbers.  (This is equivalent to the COMPLEX*16 data
*     type that most FORTRAN compilers provide but which is not part
*     of the standard language.)

*     Input:   AI, AJ, AN       The first matrix.
*              BI, BJ, BN       The second matrix.
*              CI, CJ           The structure of the result.
*              N                The number of rows.
*              M                The number of columns.

*     Output:  CN               The resulting matrix.
*              W                Workspace of length M.

      INTEGER I, J, K, CI1, CI2

      DO 50 I = 1, N

         CI1 = CI(I)
         CI2 = CI(I+1) - 1

         IF ( CI1 .GT. CI2 ) GO TO 50

         DO 10 K = CI1, CI2
            W(2*CJ(K)-1) = 0.0D0
            W(2*CJ(K)) = 0.0D0
 10      CONTINUE

         DO 20 K = AI(I), AI(I+1) - 1
            W(2*AJ(K)-1) = AN(2*K-1)
            W(2*AJ(K)) = AN(2*K)
 20      CONTINUE

         DO 30 K = BI(I), BI(I+1) - 1
            J = BJ(K)
            W(2*J-1) = W(2*J-1) + BN(2*K-1)
            W(2*J) = W(2*J) + BN(2*K)
 30      CONTINUE

         DO 40 K = CI1, CI2
            CN(2*K-1) = W(2*CJ(K)-1)
            CN(2*K) = W(2*CJ(K))
 40      CONTINUE

 50   CONTINUE

      RETURN
      END

***********************************************************************

      SUBROUTINE ZGSMUL( AI, AJ, AN, BI, BJ, BN, NA, MB,
     $     CI, CJ, CN, W )
      INTEGER AI(*), AJ(*), BI(*), BJ(*), NA, MB
      INTEGER CI(*), CJ(*)
      DOUBLE PRECISION AN(*), BN(*), CN(*), W(MB)

*     ZGSMUL multiplies two matrices in the full, unordered, sparse
*     form.  The structure of the result must have already been formed,
*     likely by calling XGSMUL.

*     Input:   AI, AJ, AN       The first matrix.
*              BI, BJ, BN       The second matrix.
*              CI, CJ           The results matrix structure.
*              NA               The number of rows of A.
*              MB               The number of columns of B.

*     Output:  CN               The resulting matrix.
*              W                Workspace of length MB.

      INTEGER I, J, K, KK, CI1, CI2
      DOUBLE PRECISION ARE, AIM, DTMP

      DO 50 I = 1, NA

         CI1 = CI(I)
         CI2 = CI(I+1) - 1

         IF ( CI1 .GT. CI2 ) GO TO 50

         DO 10 J = CI1, CI2
            W(2*CJ(J)-1) = 0.0D0
            W(2*CJ(J)) = 0.0D0
 10      CONTINUE

         DO 30 J = AI(I), AI(I+1) - 1

            ARE = AN(2*J-1)
            AIM = AN(2*J)

            DO 20 KK = BI(AJ(J)), BI(AJ(J)+1) - 1
               K    = BJ(KK)
               DTMP     = ARE*BN(2*KK-1)-AIM*BN(2*KK)
               W(2*K-1) = W(2*K-1) + DTMP
               DTMP     = ARE*BN(2*KK)+AIM*BN(2*KK-1)
               W(2*K)   = W(2*K) + DTMP
 20         CONTINUE

 30      CONTINUE

         DO 40 J = CI1, CI2
            CN(2*J-1) = W(2*CJ(J)-1)
            CN(2*J) = W(2*CJ(J))
 40      CONTINUE

 50   CONTINUE

      RETURN
      END

***********************************************************************

      SUBROUTINE ZGSSUB( AI, AJ, AN, BI, BJ, BN, N, M,
     $     CI, CJ, CN, W )
      INTEGER AI(*), AJ(*), BI(*), BJ(*), N, M
      INTEGER CI(*), CJ(*)
      DOUBLE PRECISION AN(*), BN(*), CN(*), W(*)

*     ZGSSUB subtracts double precision complex matrices in the full,
*     unordered, sparse form.  The CI and CJ vectors should be
*     previously set by a call to XGSADD.  This routine does not change
*     the order of CI and CJ, so one could call XGSORD to order them
*     before calling this routine.

*     Input:   AI, AJ, AN       The first matrix.
*              BI, BJ, BN       The second matrix.
*              CI, CJ           The structure of the result.
*              N                The number of rows.
*              M                The number of columns.

*     Output:  CN               The resulting matrix.
*              W                Workspace of length M.

      INTEGER I, J, K, CI1, CI2

      DO 50 I = 1, N

         CI1 = CI(I)
         CI2 = CI(I+1) - 1

         IF ( CI1 .GT. CI2 ) GO TO 50

         DO 10 J = CI1, CI2
            W(2*CJ(J)-1) = 0.0D0
            W(2*CJ(J)) = 0.0D0
 10      CONTINUE

         DO 20 J = AI(I), AI(I+1) - 1
            W(2*AJ(J)-1) = AN(2*J-1)
            W(2*AJ(J)) = AN(2*J)
 20      CONTINUE

         DO 30 J = BI(I), BI(I+1) - 1
            K = BJ(J)
            W(2*K-1) = W(2*K-1) - BN(2*J-1)
            W(2*K) = W(2*K) - BN(2*J)
 30      CONTINUE

         DO 40 J = CI1, CI2
            CN(2*J-1) = W(2*CJ(J)-1)
            CN(2*J) = W(2*CJ(J))
 40      CONTINUE

 50   CONTINUE

      RETURN
      END

***********************************************************************

      SUBROUTINE ZGSTRN( AI, AJ, AN, N, M, ATI, ATJ, ATN )
      INTEGER AI(*), AJ(*), ATI(*), ATJ(*), N, M
      DOUBLE PRECISION AN(*), ATN(*)

*     DGSTRN transposes a double precision complex matrix in the
*     full, unordered sparse form.

*     Input:   AI, AJ, AN       The original matrix.
*              N                The number of rows.
*              M                The number of columns.

*     Output:  ATI, ATJ, ATN    The transposed matrix.

      INTEGER I, J, K, JP

      DO 10 I = 2, M+1
         ATI(I) = 0
 10   CONTINUE

      DO 20 I = 1, AI(N+1) - 1
         J = AJ(I) + 2
         IF ( J .LE. M+1 ) ATI(J) = ATI(J) + 1
 20   CONTINUE

      ATI(1) = 1
      ATI(2) = 1

      DO 30 I = 3, M+1
         ATI(I) = ATI(I) + ATI(I-1)
 30   CONTINUE

      DO 50 I = 1, N
         DO 40 JP = AI(I), AI(I+1) - 1
            J = AJ(JP) + 1
            K = ATI(J)
            ATJ(K) = I
            ATN(2*K-1) = AN(2*JP-1)
            ATN(2*K) = AN(2*JP)
            ATI(J) = K + 1
 40      CONTINUE
 50   CONTINUE

      RETURN
      END

***********************************************************************
