      DOUBLE PRECISION FUNCTION DLARAN( ISEED )
*
*  -- LAPACK auxiliary routine (version 3.1) --
*     Univ. of Tennessee, Univ. of California Berkeley and NAG Ltd..
*     November 2006
*  -- RLABPLUS modification, September 2007
*
*     .. Array Arguments ..
      INTEGER            ISEED( 4 )

      DOUBLE PRECISION   GSLRNGUF
      EXTERNAL           GSLRNGUF

*     ..
*
*  Purpose
*  =======
*
*  GSLRNGUF returns a random real number in [0,1] from dedicated GSL RNG
*

      DLARAN = GSLRNGUF()
      RETURN
*
*     End of DLARAN
*
      END
