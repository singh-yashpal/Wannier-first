      SUBROUTINE BANDW
      INCLUDE 'PARAMS'
      INCLUDE 'commons.inc'

      COMMON/PTRANS/TVEC(3,3),atheta,iraxis,
     &               NXYZ,ISHELLX,ISHELLY,ISHELLZ
      DIMENSION XKVEC(1000), XK(1000), YK(1000), ZK(1000), TV(3)
      DIMENSION RWORK(300), EVALS(100),WORKK(1000)
      COMPLEX*16, ALLOCATABLE ::  HAMC(:,:), OVEC(:,:), WORK(:)

! === Dynamically allocated real matrices ===
      REAL*8, ALLOCATABLE :: OVERR(:,:),OVERI(:,:),HAMRE(:,:),HAMIM(:,:)
      REAL*8, ALLOCATABLE :: HPLUS(:,:,:),SPLUS(:,:,:)

      PI = 4.0D0 * DATAN(1.0D0)

      call system('rm KINBABY')
      call system('rm OVLBABY')
      call system('rm HAMBABY')
      call system('rm HNLBABY')
      call system('rm HAMOLD')

      print *,'in bandw'
      NCALC = 1
      IERR = 0
      CALL SYMBOL(NCALC,MODE_RUN,IMESH,IERR)
      OPEN(7,FILE='OUTPUT',FORM='FORMATTED',STATUS='UNKNOWN')
      REWIND(7)
      CALL FGMAT
      CALL CREPMAT
      CALL READINP


      call reset_nbas

      NBAS = NS_TOT(1)

      call nbasiscc(NCENT)
   
      PRINT *, 'NBAS IN BANDW:', NBAS
      print *,'NCENT',NCENT


      OPEN(12,FILE='KPATH.DAT',STATUS='OLD')
      READ(12,*) NK,NBAND
      READ(12,*)
      PRINT *,'NUMBER OF KPOINTS',NK
      DO IK = 1, NK
         READ(12,*) XK(IK), YK(IK), ZK(IK)
         PRINT *, XK(IK), YK(IK), ZK(IK)
      END DO
      CLOSE(12)

      ALLOCATE(HPLUS(NBAS,NBAS,1000), SPLUS(NBAS,NBAS,1000))
      ALLOCATE(OVERR(NBAS,NBAS), HAMRE(NBAS,NBAS))
      ALLOCATE(OVERI(NBAS,NBAS), HAMIM(NBAS,NBAS))
      ALLOCATE(HAMC(NBAS,NBAS), OVEC(NBAS,NBAS))
      LWORK=3*NBAS
      ALLOCATE(WORK(LWORK))

      ISTSCF=2
      DO ISPNX = 1,1!NSPN
         IF (ISPNX .EQ. 1) THEN
         OPEN(31, FILE='BANDS_UP', STATUS='UNKNOWN', FORM='FORMATTED')
         ELSE
         OPEN(31, FILE='BANDS_DN', STATUS='UNKNOWN', FORM='FORMATTED')
         ENDIF
         WRITE(31,'(A)') '# NK NBANDS NSPIN'
         WRITE(31,'(I5,1X,I5,1X,I5)') NK, NBAND, ISPNX

         HPLUS(:,:,:)  = 0.0D0
         SPLUS(:,:,:)  = 0.0D0
         HAMRE(:,:)=0.0d0
         OVERR(:,:)=0.0d0

         TV(:)     = 0.0D0
         RANG      = 0.0D0
         HSTOR(:,:)= 0.0D0

         CALL OVERFULL(1,TV,RANG)
         DO ISPN = 1, NSPN
         CALL OVERFULL(2,TV,RANG)
         END DO

         IREC = 0
         DO JBAS = 1, NBAS
         DO IBAS = 1, NBAS
            IREC = IREC + 1
            HPLUS(IBAS,JBAS,1) = HSTOR(IREC,2)
            SPLUS(IBAS,JBAS,1) = HSTOR(IREC,1)
         END DO
         END DO

         do ix=1,NBAS
         do iy=1,NBAS
         HAMRE(iy,ix)=HPLUS(iy,ix,1)
         OVERR(iy,ix)=SPLUS(iy,ix,1)
         enddo
         enddo

!         print *,'HAM'
!         do ix=1,NBAS
!         print '(14f12.5)',(HAMRE(iy,ix),iy=1,NBAS)
!         enddo
!
!         print *,'OVERLAP'
!         do ix=1,NBAS
!         print '(14f12.5)',(OVERR(iy,ix),iy=1,NBAS)
!         enddo
!
!         MDH=4*NCENT
!         CALL DSYGV(1,'N','L',NCENT,HAMRE,NBAS,OVERR,NBAS,
!     &      EVALS,WORKK,MDH,INFO)
!
!          PRINT *, 'DSYGV INFO = ', INFO
!          IF (INFO.NE.0) THEN
!             PRINT *, 'DSYGV FAILED'
!             CALL STOPIT
!          END IF
!          PRINT *, 'EVALS',(EVALS(jk),jk=1,NCENT)
!
!         call stopit

         OPEN(20, FILE='ZTRANS', STATUS='OLD')
         REWIND(20)

         DO irt=2,10000

            READ(20,*,END=1011)nx,ny,nz
            DO j=1,3
               TV(j)=nx*TVEC(j,1)+ny*TVEC(j,2)+nz*TVEC(j,3)
            ENDDO

            if(iraxis.eq.1)RANG=nx*atheta
            if(iraxis.eq.2)RANG=ny*atheta
            if(iraxis.eq.3)RANG=nz*atheta

            HSTOR(:,:) = 0.0D0
            CALL OVERFULL(1,TV,RANG)
            DO ISPN = 1, NSPN
            CALL OVERFULL(2,TV,RANG)
            END DO

            IREC = 0
            DO JBAS = 1, NBAS
            DO IBAS = 1, NBAS
               IREC = IREC + 1
               HPLUS(IBAS,JBAS,irt) = HSTOR(IREC,2)
               SPLUS(IBAS,JBAS,irt) = HSTOR(IREC,1)
            END DO
            END DO

         END DO !ir UTRANS loop
 1011    CONTINUE

!=========================================================

         dist=0.0d0
         DO IK=1,NK

            IF(IK.GT.1)THEN
            delk=dsqrt((XK(IK)-XK(IK-1))**2
     &                +(YK(IK)-YK(IK-1))**2
     &                +(ZK(IK)-ZK(IK-1))**2)
            ELSE
            delk=dsqrt(XK(IK)**2+YK(IK)**2+ZK(IK)**2)
            ENDIF

            HAMRE(:,:) = 0.0D0
            HAMIM(:,:) = 0.0D0
            OVERR(:,:) = 0.0D0
            OVERI(:,:) = 0.0D0

            DO JBAS = 1, NCENT
            DO IBAS = 1, NCENT
               HAMRE(IBAS,JBAS) = HPLUS(IBAS,JBAS,1)
               OVERR(IBAS,JBAS) = SPLUS(IBAS,JBAS,1)
            END DO
            END DO


            REWIND(20)

            DO irt=2,10000
               READ(20,*,END=1012)nx,ny,nz

               DO j=1,3
               TV(j)=nx*TVEC(j,1)+ny*TVEC(j,2)+nz*TVEC(j,3)
               ENDDO

               DOTK_R = XK(IK)*TV(1) + YK(IK)*TV(2) + ZK(IK)*TV(3)

               COSKR = DCOS(DOTK_R)
               SINKR = DSIN(DOTK_R)

               DO JBAS = 1, NCENT
               DO IBAS = 1, NCENT

                  HAMRE(IBAS,JBAS) = HAMRE(IBAS,JBAS)+
     &            COSKR*(HPLUS(IBAS,JBAS,irt))
            
                  HAMIM(IBAS,JBAS) = HAMIM(IBAS,JBAS)+
     &            SINKR*(HPLUS(IBAS,JBAS,irt))
            
                  OVERR(IBAS,JBAS) = OVERR(IBAS,JBAS)+
     &            COSKR*(SPLUS(IBAS,JBAS,irt))
            
                  OVERI(IBAS,JBAS) = OVERI(IBAS,JBAS)+
     &            SINKR*(SPLUS(IBAS,JBAS,irt))

               END DO
               END DO

            END DO !irt

 1012       continue

!
            DO JBAS = 1, NCENT
            DO IBAS = 1, NCENT
            HAMC(IBAS,JBAS) = DCMPLX(HAMRE(IBAS,JBAS),HAMIM(IBAS,JBAS))
            OVEC(IBAS,JBAS) = DCMPLX(OVERR(IBAS,JBAS),OVERI(IBAS,JBAS))
            END DO
            END DO

             ERRS = 0.0D0
             DO I=1,NCENT
             DO J=1,NCENT
                ERRS = MAX(ERRS, ABS(OVEC(I,J)-DCONJG(OVEC(J,I))))
             END DO
             END DO
!             PRINT *, 'MAX HERMITICITY ERROR S = ', ERRS

! === Enforce Hermitian symmetry on overlap matrix ===
            DO IBAS = 1, NCENT
            DO JBAS = IBAS + 1, NCENT
              OVEC(JBAS, IBAS) = DCONJG(OVEC(IBAS, JBAS))
            END DO
            END DO
! === Enforce Hermitian symmetry on Hamiltonian matrix ===
            DO IBAS = 1, NCENT
            DO JBAS = IBAS + 1, NCENT
              HAMC(JBAS, IBAS) = DCONJG(HAMC(IBAS, JBAS))
            END DO
            END DO 

            !CALL DSYGV(1,'V','L',MBAS,HAMW,NDH,VER,NDH,EVAL,SC1,MDH,info)

!            MDH=4*NCENT
!            CALL DSYGV(1,'N','U',NCENT,HAMRE,NCENT,OVERR,NCENT,
!     &          EVALS,WORKK,MDH,INFO)
            CALL ZHEGV(1,'N','U',NCENT,HAMC,NBAS,
     &               OVEC,NBAS,EVALS,WORK,LWORK,RWORK,INFO)


            IF (INFO .NE. 0) THEN
            PRINT *, 'ZHEGV failed at k =', IK, ' INFO =', INFO
            CALL STOPIT
            END IF

            !dist=(IK-1)*delk
            dist=dist+delk
            WRITE(31,'(15F15.5)') dist, (EVALS(I), I=1, 10)

         END DO !ik

!=========================================================
         CLOSE(31)
         CLOSE(20)

      END DO
      DEALLOCATE(HAMIM,HAMRE,OVERI,OVERR,HPLUS,SPLUS)
      DEALLOCATE(HAMC,OVEC,WORK)

      print *,'END OF BANDSTRUCTURE CALCULATION!'
      CALL STOPIT
      RETURN
      END
