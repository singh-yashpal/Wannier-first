C
C GET NUMERICAL PART OF HAMILTONIAN INTEGRALS
C
        IF (MODE.EQ.2)THEN
            CALL NUMHAMW(TV,RANG)
c            CALL SICSET(3)
       END IF
       RETURN
      END
C
C *****************************************************************
C
C
       SUBROUTINE NUMHAMW(TV,RANG)
C
C WRITTEN BY MARK R PEDERSON (1988-1995)
C
       INCLUDE 'PARAMS'
       INCLUDE 'commons.inc'
       LOGICAL EXIST
       COMMON/MIXPOT/POTIN(MAX_PTS*MXSPN),POTDV(MAX_PTS*MXSPN)
       DIMENSION TV(3),POTOV(MAX_PTS*MXSPN)
C
C IMODE=1 USE O(N)+PARALLEL APPROACH
C IMODE=2 USE STANDARD APPROACH - FASTER FOR SMALL SYSTEMS / FEW CPUS
C
       IMODE=2
C
C ZERO POTOLD. READ POTOLD IF POTENTIAL IS USED AS STARTING POINT
C
       CALL GTTIME(TIME1)
       DO 10 IPTS=1,NMSH*NSPN
        POTDV(IPTS)=0.0D0
        POTOV(IPTS)=1.0D0
   10  CONTINUE
       IF (ISTSCF.EQ.2) THEN
        INQUIRE(FILE='POTOLD',EXIST=EXIST)
        IF (.NOT.EXIST) THEN
         PRINT *,'NUMHAM: POTOLD DOES NOT EXIST'
         CALL STOPIT
        END IF
        OPEN(99,FILE='POTOLD',FORM='UNFORMATTED',STATUS='OLD')
        REWIND(99)
        READ(99,END=20)MPTS,MSPN,NDUM

!        IF (MPTS.NE.NMSH) THEN
!         PRINT *,'NUMHAM: NMSH IN POTOLD DIFFERS FROM CURRENT JOB'
!         PRINT*,'MPTS,NMSH=',MPTS,NMSH,MSPN,NDUM
!         CALL STOPIT
!        END IF

        NMSH=MAX(NMSH,MPTS)
        KPTS=NMSH*MIN(MSPN,NSPN)
        READ(99,END=20)(POTDV(IPTS), IPTS=1,KPTS)
        GOTO 30
   20   PRINT *,'NUMHAM: POTOLD UNREADABLE'
        CALL STOPIT
   30   CLOSE(99)

        !print *,'potpot: KPTS,NMSH',KPTS,NMSH
        !print *,(POTDV(ix),ix=1,10)
C
C DEFINE SPIN DOWN POTENTIAL IF POTOLD POTENTIAL IS SPIN UNPOLARIZED
C
        IF ((MSPN.EQ.1).AND.(NSPN.EQ.2)) THEN
         DO 40 IPTS=1,NMSH
          POTDV(IPTS+NMSH)=POTDV(IPTS)
          POTOV(IPTS+NMSH)=POTOV(IPTS)
   40    CONTINUE
        END IF

C
C MULTIPLY POTENTIAL BY WMSH AND CALL CORRECT SUBROUTINE
C
        IPTR=(ISPN-1)*NMSH

        DO IPTS=1,NMSH
         POTDV(IPTS)=WMSH(IPTS)*POTDV(IPTS+IPTR)
         POTOV(IPTS)=WMSH(IPTS)*POTOV(IPTS+IPTR)
        END DO
        IF (IMODE.EQ.1) THEN
         CALL OVERNUM(1)
        ELSE
         CALL PATCHW(POTDV,POTOV,TV,RANG)
        END IF
       END IF
C
       RETURN
       END
C
C ********************************************************************
C
C CONSTRUCT MATRIX ELEMENTS OF POTENTIAL
C MARK PEDERSON AUG 1995
C
       SUBROUTINE PATCHW(POTDV,POTOV,TV,RANG)
C
C WRITTEN BY MARK R PEDERSON
C 02/12/97 David Clay Patton
C
       INCLUDE 'PARAMS'
       INCLUDE 'commons.inc'
       PARAMETER (NMAX=MPBLOCK)
       DIMENSION POTDV(*),TV(3),POTOV(*)
C
C SCRATCH COMMON BLOCK FOR LOCAL ARRAYS
C
       COMMON/TMP2/PSIBR(MPBLOCK,5*NDH),PSIB(MPBLOCK)
     &  ,VOL(MPBLOCK),QR(3,MPBLOCK)
       common/tmp22/PSIBL(MPBLOCK,5*NDH),TR(3,MPBLOCK)

       LOGICAL LSETUP
       DATA MDCL/2/
C
C GAM=   REPRESENATION INDEX, L=ROW OF REPRESENATION
C MDCL=1 CALCULATE MATRIX ELEMENTS BY ROTATING EACH MESH POINT TO
C        EQUIVALENT POINTS AND FIND <PSI(I,GAM,L)|V(I)|PSI(I,GAM,L)>
C        SUMMED OVER I WITH L=1
C MDCL=2 CALCULATE MATRIX ELEMENTS BY USING ONLY INEQUVIALENT POINTS
C        AND FIND <PSI(I,GAM,L)|V(I)|PSI(I,GAM,L)> SUMMED OVER I AND L
C
       IF (NMSH.GT.MAX_PTS) THEN
        PRINT *,'PATCH: MAX_PTS MUST BE AT LEAST', NMSH
        CALL STOPIT
       END IF
C
C      INITIALIZE ARRAY INDEX IN GETBAS
C
       LSETUP=.TRUE.
       CALL GETBAS(LSETUP,0,0,QR,PSIBR,0,NBAS)
       CALL GETBAS(LSETUP,0,0,TR,PSIBL,0,NBAS)
C
       DO LPTS=0,NMSH-1,NMAX
        MPTS=MIN(NMAX,NMSH-LPTS)
        CALL HAMMATW(LPTS,MPTS,MDCL,POTDV,POTOV,INDX,TV,RANG)
       END DO
       RETURN
       END
C
      SUBROUTINE HAMMATW(LPTS,MPTS,MDCL,POTDV,POTOV,INDX,TV,RANG)
C
C 02/12/97 David Clay Patton
C
      INCLUDE 'PARAMS'
      INCLUDE 'commons.inc'
      PARAMETER (pi = 3.14159265358979323846d0)
      LOGICAL LSETUP,LSETDN
      DIMENSION POTDV(*),TV(3),POTOV(*),VOLOV(MPBLOCK)
      DIMENSION vin(3),vout(3),PSILAP(MPBLOCK,5*NDH)
      DIMENSION PSIKE(MPBLOCK)
      COMMON/TMP2/PSIBR(MPBLOCK,5*NDH),PSIB(MPBLOCK),PSIBOV(MPBLOCK)
     & ,VOL(MPBLOCK),QR(3,MPBLOCK)
      common/tmp22/PSIBL(MPBLOCK,5*NDH),TR(3,MPBLOCK)
      COMMON/PTRANS/TVEC(3,3),atheta,iraxis,
     &               NXYZ,ISHELLX,ISHELLY,ISHELLZ      
     
C

      MGRP=1
      IF (MDCL.EQ.1) MGRP=NGRP
C
      LSETUP=.FALSE.
      DO IPTS=1,MPTS

         VOL(IPTS)=POTDV(LPTS+IPTS)/MGRP
         VOLOV(IPTS)=POTOV(LPTS+IPTS)/MGRP

         QR(1,IPTS)=RMSH(1,LPTS+IPTS)
         QR(2,IPTS)=RMSH(2,LPTS+IPTS)
         QR(3,IPTS)=RMSH(3,LPTS+IPTS)

         vin(1)=QR(1,IPTS)
         vin(2)=QR(2,IPTS)
         vin(3)=QR(3,IPTS)

         if(abs(RANG).gt.0.d-10)then
         call rotvec(iraxis,RANG,vin,vout)
         else
         vout(1)=vin(1)
         vout(2)=vin(2)
         vout(3)=vin(3)
         endif

         TR(1,IPTS)=vout(1)+TV(1)
         TR(2,IPTS)=vout(2)+TV(2)
         TR(3,IPTS)=vout(3)+TV(3)

      END DO
C
C UPDATE MATRIX ELEMENTS HERE:
C
      DO 100 IGP=1,MGRP
       JGP=IGP
       INDX=0
       DO 200 KREP=1,N_REP
        JREP=KREP
        IF (MDCL.EQ.1) THEN
         NRP=1
        ELSE
         NRP=NDMREP(KREP)
        END IF

        PSILAP = 0.0d0
        CALL GETBAS(LSETUP,MPTS,JGP,QR,PSIBR,JREP,NBAS)
        CALL GETBAS(LSETDN,MPTS,JGP,TR,PSIBL,JREP,NBAS)
        !CALL GETLAP(MPTS,JGP,QR,PSILAP,JREP,NBAS)

        JNDX=INDX
        DO 300 KRP=1,NRP
         INDX=JNDX
         IBG=(KRP-1)*NBAS
         DO 400 ISS=1,NBAS
          DO IPTS=1,MPTS
           PSIB(IPTS)  = PSIBR(IPTS,ISS+IBG) *  VOL(IPTS)/NRP
           PSIBOV(IPTS)= PSIBR(IPTS,ISS+IBG) *VOLOV(IPTS)/NRP
           PSIKE(IPTS) = PSILAP(IPTS,ISS+IBG)*VOLOV(IPTS)/NRP
          END DO
C NOTE: WE COULD MAKE USE OF SPARCITY HERE. IF PSIB=0 AT EACH
C       MESHPOINT, THE NEXT LOOP CAN BE SKIPPED - INDX MUST BE
C       RESET OF COURSE
C
          DO JSS=ISS,NBAS
             DOT=0.0D0
             DOTOV=0.0D0
             DOTKE=0.0D0
             DO IPTS=1,MPTS
                DOT  =DOT+  PSIB(IPTS)*PSIBL(IPTS,JSS+IBG)
                DOTOV=DOTOV+PSIBOV(IPTS)*PSIBL(IPTS,JSS+IBG)
                DOTKE=DOTKE+PSIKE(IPTS)*PSIBL(IPTS,JSS+IBG)
             END DO
             INDX=INDX+1
             HSTOR(INDX,2)=HSTOR(INDX,2)+(DOT+DOTKE)
             !HSTOR(INDX,1)=HSTOR(INDX,1)+ DOTOV
          END DO
400      CONTINUE
300     CONTINUE
200    CONTINUE
100   CONTINUE

      RETURN
      END
!========================================================

      SUBROUTINE GETLAP(MPTS,JGP,QR,PSITAP,JREP,NBAS)
       INCLUDE 'PARAMS'
       INCLUDE 'commons.inc'
       PARAMETER (NMAX=MPBLOCK)
C
       LOGICAL ICOUNT
       COMMON/TMP2/PSIG(NMAX,10,MAX_OCC)
     &  ,PTS(NSPEED,3),GRAD(NSPEED,10,6,MAX_CON,3)
     &  ,RVECA(3,MX_GRP),ICOUNT(MAX_CON,3)
C
C SCRATCH COMMON BLOCK FOR LOCAL ARRAYS
C
       LOGICAL LGGA,IUPDAT
       DIMENSION ISIZE(3),QR(3,MPBLOCK)
       DIMENSION PSITAP(MPBLOCK,5*NDH)
       DATA ISIZE/1,3,6/

                    ! Laplacian = d²/dx² + d²/dy² + d²/dz²
                    ! GRAD indices: (point, ?, orbital, contraction, L)
                    ! GRAD(LPV,5,MU,ICON,LI) = d²/dx²
                    ! GRAD(LPV,6,MU,ICON,LI) = d²/dy²
                    ! GRAD(LPV,7,MU,ICON,LI) = d²/dz²
                   
C
C FOR ALL CENTER TYPES
C
        IBAS=0
        PSITAP(:,:) = 0.0d0

        ISHELLA=0

        !print *,'NFNCT',NFNCT

        DO 86 IFNCT=1,NFNCT
           LMAX1=LSYMMAX(IFNCT)+1
!           PRINT *,'IFNCT=',IFNCT,' N_POS=',N_POS(IFNCT),
!     &             ' LMAX1=',LMAX1
C
C FOR ALL POSITIONS OF THIS CENTER
C
          DO 84 I_POS=1,N_POS(IFNCT)
             ISHELLA=ISHELLA+1
C
C GET SYMMETRY INFO
C
             CALL OBINFO(1,RIDT(1,ISHELLA),RVECA,M_NUC,ISHDUM)
!             PRINT *,'  I_POS=',I_POS,' M_NUC=',M_NUC
C
C FOR ALL EQUIVALENT POSITIONS OF THIS ATOM
C
            DO 82 J_POS=1,M_NUC
C
C UNSYMMETRIZE 
C
            CALL UNRAVEL(IFNCT,ISHELLA,J_POS,RIDT(1,ISHELLA),
     &                  RVECA,L_NUC,1)
            IF(L_NUC.NE.M_NUC)THEN
               PRINT *,'APTSLV: PROBLEM IN UNRAVEL'
            CALL STOPIT
            END IF
C
C FOR ALL MESHPOINTS IN BLOCK DO A SMALLER BLOCK
C
            IBAS_START=IBAS

            KPTS=0
            DO 80 JPTS=1,MPTS,NSPEED
               NPV=MIN(NSPEED,MPTS-JPTS+1)
               DO LPV=1,NPV
                  KPTS=KPTS+1
                  PTS(LPV,1)=QR(1,KPTS)-RVECA(1,J_POS)
                  PTS(LPV,2)=QR(2,KPTS)-RVECA(2,J_POS)
                  PTS(LPV,3)=QR(3,KPTS)-RVECA(3,J_POS)
               END DO
C
C GET ORBITS AND DERIVATIVES
C
              NDERV=2
              CALL GORBDRV(NDERV,IUPDAT,ICOUNT,NPV,PTS,IFNCT,GRAD)

C
C ACCUMULATE LAPLACIAN: -0.5 * ∇²ψ
C
              IF(IUPDAT)THEN
              IPTS=JPTS-1
              IBAS=IBAS_START
              DO LI=1,LMAX1
              ! ICON loop before MU: 2px 2py 2pz, 3px 3py, 3pz  and so on 
              ! MU loop before ICON 2px 3px 4px, 2py 3py 4py and so on
                 DO ICON=1,N_CON(LI,IFNCT)
                 DO MU=1,ISIZE(LI)
                    IBAS=IBAS+1
                    IF(ICOUNT(ICON,LI))THEN
                    DO LPV=1,NPV
                       PSITAP(IPTS+LPV,IBAS)= -0.5D0*
     &                        (GRAD(LPV,5,MU,ICON,LI)
     &                   +     GRAD(LPV,6,MU,ICON,LI)
     &                   +     GRAD(LPV,7,MU,ICON,LI))
                    END DO
                    ENDIF
                 END DO
                 END DO
              END DO

           ENDIF
   80      CONTINUE
           IBAS = IBAS_START
           DO LI=1,LMAX1
              IBAS = IBAS + ISIZE(LI)*N_CON(LI,IFNCT)
           END DO
!           PRINT *,'    J_POS=',J_POS,' basis range:',
!     &     IBAS_START+1,' to ',IBAS

   82     CONTINUE
   84     CONTINUE
   86     CONTINUE

        RETURN
        END

