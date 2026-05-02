      SUBROUTINE BANDWMPI
      INCLUDE 'PARAMS'
      INCLUDE 'commons.inc'
      INCLUDE 'mpif.h'

      COMMON /PTRANS/ TVEC(3,3), ATHETA, IPAX, NXYZ, ISHELLX
      DIMENSION XKVEC(1000), XK(100), YK(100), ZK(100), TV(3)
      DIMENSION RWORK(300), EVALS(100),WORKK(1000)
      COMPLEX*16, ALLOCATABLE ::  HAMC(:,:), OVEC(:,:), WORK(:)

! === Dynamically allocated real matrices ===
      REAL*8, ALLOCATABLE :: OVERR(:,:),OVERI(:,:),HAMRE(:,:),HAMIM(:,:)
      REAL*8, ALLOCATABLE :: HPLUS(:,:,:),SPLUS(:,:,:)

! === MPI variables ===
      INTEGER :: MYRANK, NPROCS, IERR, NELEM

! === Pre-read arrays for UTRANS ===
      INTEGER, ALLOCATABLE :: NXA(:), NYA(:), NZA(:)
      INTEGER :: NTRANS, NTOT

      PI = 4.0D0 * DATAN(1.0D0)
      NBAS = NS_TOT(1)
      PRINT *, 'NBAS IN BANDW:', NBAS
      call nbasiscc(NCENT)

      print *,'NCENT NSPN',NCENT,NPSN

      CALL MPI_COMM_RANK(MPI_COMM_WORLD, MYRANK, IERR)
      CALL MPI_COMM_SIZE(MPI_COMM_WORLD, NPROCS, IERR)

      OPEN(12,FILE='KPATH.DAT',STATUS='OLD')
      READ(12,*) NK,NBAND
      READ(12,*)
      IF (MYRANK .EQ. 0) PRINT *,'NUMBER OF KPOINTS',NK
      DO IK = 1, NK
         READ(12,*) XK(IK), YK(IK), ZK(IK)
         IF (MYRANK .EQ. 0) PRINT *, XK(IK), YK(IK), ZK(IK)
      END DO
      CLOSE(12)

C ============================================================
C  Pre-read ALL translations from UTRANS into arrays
C  (every rank reads — small file, avoids broadcast)
C ============================================================
      ALLOCATE(NXA(10000), NYA(10000), NZA(10000))
      OPEN(20, FILE='UTRANS', STATUS='OLD')
      NTRANS = 0
      DO irt = 1, 10000
         READ(20,*,END=900) NXA(irt), NYA(irt), NZA(irt)
         NTRANS = irt
      END DO
  900 CONTINUE
      CLOSE(20)
      IF (MYRANK .EQ. 0) PRINT *, 'Read NTRANS translations:', NTRANS

C --- Total slots: home cell (1) + NTRANS translations ---
      NTOT = NTRANS + 1

      ALLOCATE(HPLUS(NCENT,NCENT,NTOT), SPLUS(NCENT,NCENT,NTOT))
      ALLOCATE(OVERR(NCENT,NCENT), HAMRE(NCENT,NCENT))
      ALLOCATE(OVERI(NCENT,NCENT), HAMIM(NCENT,NCENT))
      ALLOCATE(HAMC(NCENT,NCENT), OVEC(NCENT,NCENT))
      LWORK=3*NCENT
      ALLOCATE(WORK(LWORK))

      print *,NSPN
      ISTSCF=2
      DO ISPNX = 1,1

         IF (MYRANK .EQ. 0) THEN
            IF (ISPNX .EQ. 1) THEN
            OPEN(31,FILE='BANDS_UP',STATUS='UNKNOWN',FORM='FORMATTED')
            ELSE
            OPEN(31,FILE='BANDS_DN',STATUS='UNKNOWN',FORM='FORMATTED')
            ENDIF
            WRITE(31,'(A)') '# NK NBANDS NSPIN'
            WRITE(31,'(I5,1X,I5,1X,I5)') NK, NBAND, ISPNX
         ENDIF

         HPLUS(:,:,:)  = 0.0D0
         SPLUS(:,:,:)  = 0.0D0
         HAMRE(:,:)    = 0.0D0
         OVERR(:,:)    = 0.0D0

C ============================================================
C  Home cell (slot 1): TV=0, RANG=0
C  Every rank computes this (cheap, avoids broadcast)
C ============================================================
         TV(:)     = 0.0D0
         RANG      = 0.0D0
         HSTOR(:,:)= 0.0D0

         CALL OVERLAPT(1,TV,RANG)
         DO ISPN = 1, NSPN
            CALL OVERLAPT(2,TV,RANG)
         END DO

         IREC = 0
         DO IBAS = 1, NCENT
         DO JBAS = IBAS, NCENT
            IREC = IREC + 1
            HPLUS(IBAS,JBAS,1) = HSTOR(IREC,2)
            SPLUS(IBAS,JBAS,1) = HSTOR(IREC,1)
         END DO
            IREC = IREC + (NBAS-NCENT)
         END DO

C ============================================================
C  ALL ranks call OVERLAPT for every irt (keeps MPI collectives
C  inside OVERLAPT synchronized). Only the owning rank stores
C  the result into HPLUS/SPLUS.
C ============================================================
         DO irt = 1, NTRANS

            nx = NXA(irt)
            ny = NYA(irt)
            nz = NZA(irt)

            DO j = 1, 3
               TV(j) = nx*TVEC(j,1) + ny*TVEC(j,2) + nz*TVEC(j,3)
            ENDDO
            RANG = nx * atheta

            HSTOR(:,:) = 0.0D0
            CALL OVERLAPT(1, TV, RANG)
            DO ISPN = 1, NSPN
               CALL OVERLAPT(2, TV, RANG)
            END DO

C --- Only the owning rank stores into HPLUS/SPLUS ---
            IF (MOD(irt-1, NPROCS) .EQ. MYRANK) THEN
               IREC = 0
               DO IBAS = 1, NCENT
               DO JBAS = IBAS, NCENT
                  IREC = IREC + 1
                  HPLUS(IBAS,JBAS,irt+1) = HSTOR(IREC,2)
                  SPLUS(IBAS,JBAS,irt+1) = HSTOR(IREC,1)
               END DO
                  IREC = IREC + (NBAS-NCENT)
               END DO
            ENDIF

         END DO

C --- Gather all pieces ---
         NELEM = NCENT * NCENT * NTOT

         CALL MPI_ALLREDUCE(MPI_IN_PLACE, HPLUS, NELEM,
     &        MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, IERR)

         CALL MPI_ALLREDUCE(MPI_IN_PLACE, SPLUS, NELEM,
     &        MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, IERR)


         IF (MYRANK .EQ. 0) THEN
            PRINT *, 'HPLUS/SPLUS built and synchronized across ranks.'
         ENDIF

C =========================================================
C  K-point loop: all ranks compute, only rank 0 writes
C =========================================================
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

C --- Home cell contribution ---
            DO IBAS = 1, NCENT
            DO JBAS = IBAS, NCENT
               HAMRE(IBAS,JBAS) = HPLUS(IBAS,JBAS,1)
               OVERR(IBAS,JBAS) = SPLUS(IBAS,JBAS,1)
            END DO
            END DO

C --- Fourier sum over all translations ---
            DO irt = 1, NTRANS

               DO j = 1, 3
               TV(j) = NXA(irt)*TVEC(j,1) + NYA(irt)*TVEC(j,2)
     &               + NZA(irt)*TVEC(j,3)
               ENDDO

               DOTK_R = XK(IK)*TV(1) + YK(IK)*TV(2) + ZK(IK)*TV(3)
               COSKR = DCOS(DOTK_R)
               SINKR = DSIN(DOTK_R)

               DO IBAS = 1, NCENT
               DO JBAS = IBAS, NCENT

                  HAMRE(IBAS,JBAS) = HAMRE(IBAS,JBAS)+
     &            COSKR*(HPLUS(IBAS,JBAS,irt+1))

                  HAMIM(IBAS,JBAS) = HAMIM(IBAS,JBAS)+
     &            SINKR*(HPLUS(IBAS,JBAS,irt+1))

                  OVERR(IBAS,JBAS) = OVERR(IBAS,JBAS)+
     &            COSKR*(SPLUS(IBAS,JBAS,irt+1))

                  OVERI(IBAS,JBAS) = OVERI(IBAS,JBAS)+
     &            SINKR*(SPLUS(IBAS,JBAS,irt+1))

               END DO
               END DO

            END DO

C --- Build complex matrices ---
            DO IBAS = 1, NCENT
            DO JBAS = IBAS, NCENT
            HAMC(IBAS,JBAS)=DCMPLX(HAMRE(IBAS,JBAS),HAMIM(IBAS,JBAS))
            OVEC(IBAS,JBAS)=DCMPLX(OVERR(IBAS,JBAS),OVERI(IBAS,JBAS))
            END DO
            END DO

C === Enforce Hermitian symmetry on overlap matrix ===
            DO IBAS = 1, NCENT
            DO JBAS = IBAS + 1, NCENT
              OVEC(JBAS, IBAS) = DCONJG(OVEC(IBAS, JBAS))
            END DO
            END DO

C === Enforce Hermitian symmetry on Hamiltonian matrix ===
            DO IBAS = 1, NCENT
            DO JBAS = IBAS + 1, NCENT
              HAMC(JBAS, IBAS) = DCONJG(HAMC(IBAS, JBAS))
            END DO
            END DO

C --- Diagonalize ---
            CALL ZHEGV(1,'N','U',NCENT,HAMC,NCENT,
     &             OVEC,NCENT,EVALS,WORK, LWORK,RWORK,INFO)

            IF (INFO .NE. 0) THEN
            PRINT *, 'ZHEGV failed at k =', IK, ' INFO =', INFO
            CALL STOPIT
            END IF

            dist=dist+delk

C --- Only rank 0 writes output ---
            IF (MYRANK .EQ. 0) THEN
               WRITE(31,'(15F15.5)') dist, (EVALS(I), I=1, 10)
            ENDIF

         END DO !ik

C =========================================================
         IF (MYRANK .EQ. 0) CLOSE(31)

      END DO

      DEALLOCATE(HAMIM,HAMRE,OVERI,OVERR,HPLUS,SPLUS)
      DEALLOCATE(HAMC,OVEC,WORK)
      DEALLOCATE(NXA,NYA,NZA)

      IF (MYRANK .EQ. 0) THEN
         print *,'END OF BANDSTRUCTURE CALCULATION!'
      ENDIF
      CALL STOPIT
      RETURN
      END
