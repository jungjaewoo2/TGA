c
c
c    DODRL - simple routine
c
c     ===============================================================
      SUBROUTINE DODRL
     +  (ODRFNC,N,M,NP,NQ,
     +   BETA,
     +   Y,LDY,X,LDX,
     +   WE,LDWE,LD2WE,WD,LDWD,LD2WD,
     +   JOB,
     +   odrout, lodrout,
     +   WORK,LWORK,IWORK,LIWORK,
     +   INFO)
c     ===============================================================
c
c     An ODRPACK driver for rlab. See Reference manual, pp. 25
c
c     This version by Marijan Kostrun, project rlabplus, VIII-2005

      external ODRFNC, DWINF

c      LOGICAL ISODR
c      INTEGER DELTAI,EPSI,XPLUSI,FNI,SDI,VCVI,
c     +        RVARI,WSSI,WSSDEI,WSSEPI,RCONDI,ETAI,
c     +        OLMAVI,TAUI,ALPHAI,ACTRSI,PNORMI,RNORSI,PRERSI,
c     +        PARTLI,SSTOLI,TAUFCI,EPSMAI,
c     +        BETAOI,BETACI,BETASI,BETANI,SI,SSI,SSFI,QRUAXI,UI,
c     +        FSI,FJACBI,WE1I,DIFFI,
c     +        DELTSI,DELTNI,TI,TTI,OMEGAI,FJACDI,
c     +        WRK1I,WRK2I,WRK3I,WRK4I,WRK5I,WRK6I,WRK7I,
c     +        LWKMN

      integer lodrout
      character*256 odrout

      integer IPRINT,LUNERR,LUNRPT

      if(lodrout.gt.1) then
          IPRINT = 6666
          open (unit=6, FILE=odrout, STATUS='OLD' )
      else
          IPRINT = 0
      endif

      LUNERR = -1
      LUNRPT = -1

      call dodr
     +  (ODRFNC,
     +   N,M,NP,NQ,
     +   BETA,
     +   Y,LDY,X,LDX,
     +   WE,LDWE,LD2WE,WD,LDWD,LD2WD,
     +   JOB,
     +   IPRINT,LUNERR,LUNRPT,
     +   WORK,LWORK,IWORK,LIWORK,
     +   INFO)

c      CALL DWINF
c     +   (N,M,NP,NQ,LDWE,LD2WE,ISODR,
c     +   DELTAI,EPSI,XPLUSI,FNI,SDI,VCVI,
c     +   RVARI,WSSI,WSSDEI,WSSEPI,RCONDI,ETAI,
c     +   OLMAVI,TAUI,ALPHAI,ACTRSI,PNORMI,RNORSI,PRERSI,
c     +   PARTLI,SSTOLI,TAUFCI,EPSMAI,
c     +   BETAOI,BETACI,BETASI,BETANI,SI,SSI,SSFI,QRUAXI,UI,
c     +   FSI,FJACBI,WE1I,DIFFI,
c     +   DELTSI,DELTNI,TI,TTI,OMEGAI,FJACDI,
c     +   WRK1I,WRK2I,WRK3I,WRK4I,WRK5I,WRK6I,WRK7I,
c     +   LWKMN)
c
      if(lodrout.gt.1) close(unit=6)

      return
      end


c
c
c    DODRCRL - more complex routine
c
c     ===============================================================
      SUBROUTINE DODRCRL
     &  (ODRFNC,N,M,NP,NQ,
     &   BETA,
     &   Y,LDY,X,LDX,
     &   WE,LDWE,LD2WE,WD,LDWD,LD2WD,
     &   IFIXB,IFIXX,LDIFX,
     &   JOB,NDIGIT,TAUFAC,
     &   SSTOL,PARTOL,MAXIT,
c     &   IPRINT,LUNERR,LUNRPT,
     &   odrout, lodrout,
     &   STPB,STPD,LDSTPD,
     &   SCLB,SCLD,LDSCLD,
     &   WORK,LWORK,IWORK,LIWORK,
     &   INFO)
c     ===============================================================
c
c     An ODRPACK driver for rlab. See Reference manual, pp. 25
c
c     This version by Marijan Kostrun, project rlabplus, VIII-2005

      external ODRFNC

      integer lodrout
      character*256 odrout

      integer IPRINT,LUNERR,LUNRPT

      if(lodrout.gt.1) then
          IPRINT = 6666
          open (unit=6, FILE=odrout, STATUS='OLD' )
      else
          IPRINT = 0
      endif

      LUNERR = -1
      LUNRPT = -1

      call dodrc
     &   (ODRFNC,
     &   N,M,NP,NQ,
     &   BETA,
     &   Y,LDY,X,LDX,
     &   WE,LDWE,LD2WE,WD,LDWD,LD2WD,
     &   IFIXB,IFIXX,LDIFX,
     &   JOB,NDIGIT,TAUFAC,
     &   SSTOL,PARTOL,MAXIT,
     &   IPRINT,LUNERR,LUNRPT,
     &   STPB,STPD,LDSTPD,
     &   SCLB,SCLD,LDSCLD,
     &   WORK,LWORK,IWORK,LIWORK,
     &   INFO)
c
c
      if(lodrout.gt.1) close(unit=6)
c
c
      return
      end

