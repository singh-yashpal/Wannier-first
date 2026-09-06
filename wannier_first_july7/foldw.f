      SUBROUTINE FOLDW(IFOLD)
      INCLUDE 'PARAMS'
      INCLUDE 'commons.inc'
      PARAMETER (LM=06)
      PARAMETER (LMX=3*LM)
      PARAMETER (MAX_ANG=((LMX+1)*(LMX+2))/2)
      PARAMETER (MPX=MAX_ANG)
      PARAMETER (MXPR=MXPOISS)
      DIMENSION POLY(MPX,(LM+1)**2),SMULTI((LM+1)**2)
      DIMENSION vin(3),vout(3),TV(3),RMSHCC(3,MAX_PTS)
      DIMENSION COULCC(MAX_PTS),RHOCC(MAX_PTS,NVGRAD,MXSPN)
      LOGICAL FIRST,FIRSTLAT,EXIST
      COMMON/FLOMESH/RABCD(3,1000)
      COMMON/MIXPOT/POTIN(MAX_PTS*MXSPN),POT(MAX_PTS*MXSPN)
      COMMON/TMP1/COULOMB(MAX_PTS),RHOG(MAX_PTS,NVGRAD,MXSPN)
      COMMON/REPLICA/NMSHW,NDCELL,MMSH1,NATMS,
     &      ZWATM(1000),RATMS(3,1000),ITSCF
      COMMON/PTRANS/TVEC(3,3),ATHETA,IRAXIS,NXYZ,iwx,iwy,iwz
      COMMON/COULDW/NANG,ANG(3,MAX_ANG),DOMEGANG(MAX_ANG),RAD_SPH
      COMMON/TEMPW/AMSH(3,500),ATWD(3,1000)
      DATA FIRST/.TRUE./
      DATA FIRSTLAT/.TRUE./

      pi=4.0d0*atan(1.0d0)

! NMSHW: Mesh points in the central cell.
! Calling COUPOT in APOTNL to evalaute COULOMB potential on NMSH grid of wannier
! domain
      IF(IFOLD.EQ.1)THEN
!==============================================================================
! UTRANS contain all the nx, ny, and nz associated withe Wannier domain. In
! conjuction with lattice translation vectors TVEC, here we evalaute external
! potential on grid points in central cell due to nuclie in WD as POTINW  
!==============================================================================
       CALL GTTIME(TPOT1)
!================================================================================
! V_ext(eN) from Wannier Domain evaluated in Central Cell.

       DO i=1,NMSH
          RMSHCC(1,i)=RMSH(1,i)
          RMSHCC(2,i)=RMSH(2,i)
          RMSHCC(3,i)=RMSH(3,i)
       ENDDO
       RMSH(:,:)=0.0d0

       INQUIRE(FILE='ZTRANS', EXIST=EXIST)
 
       IF (EXIST) THEN
          OPEN(20, FILE='ZTRANS', STATUS='OLD')
       ELSE
          INQUIRE(FILE='UTRANS', EXIST=EXIST)
          IF (.NOT. EXIST) STOP 'Neither ZTRANS nor UTRANS found'
          OPEN(20, FILE='UTRANS', STATUS='OLD')
       END IF
 
       REWIND(20)

       DELTA=1.0D-100
       COULCC = 0.0d0
       RHOCC = 0.0d0
       POTIN=0.0D0

       do iwxx=1,10000
          COULOMB(:)  = 0.0d0
          RHOG(:,:,:) = 0.0d0

          if(iwxx.eq.1)then
            nx=0
            ny=0
            nz=0
          else
            read(20,*,end=1011)nx,ny,nz
          endif

          do ix=1,3
             TV(ix)=nx*TVEC(ix,1)+ny*TVEC(ix,2)+nz*TVEC(ix,3)
          enddo

          do IMSH = 1,NMSH

             vin(1)=RMSHCC(1,IMSH)+TV(1)
             vin(2)=RMSHCC(2,IMSH)+TV(2)
             vin(3)=RMSHCC(3,IMSH)+TV(3)

             if(iraxis.eq.1)RANG=nx*atheta
             if(iraxis.eq.2)RANG=ny*atheta
             if(iraxis.eq.3)RANG=nz*atheta

             if(abs(RANG).gt.1.0d-10)then
             call rotvec(iraxis,RANG,vin,vout)
             else
             vout(1)=vin(1)
             vout(2)=vin(2)
             vout(3)=vin(3)
             endif

             do ix=1,3
                RMSH(ix,IMSH)=vout(ix)
             enddo

             if(iwxx.eq.1)then
             do iatm=1,NATMS*NDCELL
                vx=ATWD(1,iatm)
                vy=ATWD(2,iatm)
                vz=ATWD(3,iatm)

                DIST = (vout(1)-vx)**2+(vout(2)-vy)**2+
     &              (vout(3)-vz)**2
                DIST=MAX(DIST,DELTA)
                DIST=1.0D0/DSQRT(DIST)
                POTIN(IMSH)= POTIN(IMSH)-ZWATM(iatm)*DIST !make ZATMS general
             enddo
             endif
          enddo


       CALL COUPOT1

!       print *,'MODDEN',MODDEN
!
!       IF (MODDEN .EQ. 2) THEN
!        CALL DENSOLD(TIMEGORB)
!       END IF
!
C
C UPDATE DATA IN RHOG
C
       NGRAD=1
       IF ((IGGA(1).GT.0).OR.(IGGA(2).GT.0)) NGRAD=10
         DO IPTS=1,NMSH
         if(abs(nx).le.iwx.and.abs(ny).le.iwy.and.abs(nz).le.iwz)then
         COULCC(IPTS)=COULCC(IPTS)+COULOMB(IPTS)
         endif
         DO ISPN =1,NSPN
         DO IGRAD=1,NGRAD
         RHOCC(IPTS,IGRAD,ISPN)= RHOCC(IPTS,IGRAD,ISPN)+
     &                           RHOG(IPTS,IGRAD,ISPN)
         ENDDO
         ENDDO
        ENDDO

       enddo

 1011  continue
       CALL GTTIME(TPOT2)

       PRINT *,'TIME ELAPSED IN POTINW',(TPOT2-TPOT1)

       PRINT *,'==============================================='
       PRINT *,'  FOLDING FROM WANNIER DOMAIN TO CENTRAL CELL'
       PRINT *,'==============================================='

! CENTRAL CELL CHARGE BEFORE FOLDING

        qtest2=0.0d0
        qtest=0.0d0
        DO ipts=1,NMSH
          qtest=qtest+(RHOCC(ipts,1,1)+RHOCC(ipts,1,2))*WMSH(ipts)
        ENDDO

        print *,'Total Charge:',qtest

        DO IMSH=1,NMSH
           COULCC(IMSH)=COULCC(IMSH)+POTIN(IMSH)
        ENDDO 

        ECOUL=0.0D0
        DO IPTS=1,NMSH
           ECOUL = ECOUL+0.5D0*COULCC(IPTS)*
     &             (RHOCC(IPTS,1,1)+RHOCC(IPTS,1,2))*WMSH(IPTS)
        ENDDO 

        print *,'ECOUL',ECOUL
!==========================================================================
! Coulomb potential is evaluated at each nuclear site in the Wannier domain 
! and then integrated
!==========================================================================
       RMSH =0.0d0
       do iat=1,NDCELL*NATMS
       do ix=1,3
          RMSH(ix,iat)=ATWD(ix,iat)
       enddo
       enddo

       COULOMB(:) = 0.0d0
       CALL COUPOT1

       POTIN=0.0d0
       DO iat=1,NATMS
       DO ibt=0,NDCELL-1
          pottemp   =COULOMB(iat+ibt*NATMS)
          POTIN(iat)=POTIN(iat)+pottemp
       ENDDO
       ENDDO

       ELOCALW=0.0D0
       DO j=1,NATMS
       ELOCALW=ELOCALW-0.5D0*ZWATM(j)*POTIN(j)
       ENDDO

        print *,'ELOCALW',ELOCALW
!==========================================================================
! ENNUCW
!==========================================================================
       IF(FIRST)THEN
         PRINT *,'Nuclear coordinates in Wannier Domain'
         do iat=1,NDCELL*NATMS
         print *,(ATWD(j,iat),j=1,3)
         enddo

         POTIN=0.0d0
         DO iat=1,NATMS
         DO jat=1,NATMS
         DO iwd=0,NDCELL-1
         vx=ATWD(1,iwd*natms+jat)
         vy=ATWD(2,iwd*natms+jat)
         vz=ATWD(3,iwd*natms+jat)

         DIST = (ATWD(1,iat)-vx)**2+
     &          (ATWD(2,iat)-vy)**2+
     &          (ATWD(3,iat)-vz)**2
         IF(DIST.gt.1.0D-10)THEN
         DIST=1.0D0/DSQRT(DIST)
         POTIN(iat) = POTIN(iat)+ZWATM(jat)*DIST
         ENDIF
         ENDDO
         ENDDO
         ENDDO

! EZNUCW is the energy due to interaction of central cell nuclie with the
! nuclie's in WD (E_NN belong to WD).
         ENNUCW=0.0d0
         DO j=1,NATMS
         ENNUCW=ENNUCW+0.5d0*ZWATM(j)*POTIN(j)
         ENDDO
         FIRST=.FALSE.
         print *,'rad_sph (fold)',rad_sph
       ENDIF !FIRST

        print *,'ENNUCW',ENNUCW
!==========================================================================
! Skipping SMULTI and LATTICESUM in the first iteration.
!==========================================================================
       IF(FIRSTLAT)THEN

!          DO ix=1,NANG
!             vecin(1)=ANG(1,ix)
!             vecin(2)=ANG(2,ix)
!             vecin(3)=ANG(3,ix)
!             call rotvec(3,90,vecin,vecout)
!             ANG(1,ix)=vecout(1)
!             ANG(2,ix)=vecout(2)
!             ANG(3,ix)=vecout(3)
!          ENDDO
!
!          POLYR=0.0d0
!          CALL HARMONICS(MPX,NANG,LM,ANG,POLY,NPOLY)




         CALL HARMONICS(MPX,NANG,LM,ANG,POLY,NPOLY)

!         print *,'printing POLY'
!         do ip=1,NPOLY
!         do ix=1,NANG
!           print *,(POLY(ix,ip))
!         enddo
!         enddo

!         DO IP=MMSH1+NDCELL*NATMS+1,NMSH
!            vecin(1)=RMSH(1,IP)
!            vecin(2)=RMSH(2,IP)
!            vecin(3)=RMSH(3,IP)
!            call rotvec(1,90,vecin,vecout)
!            RMSH(1,IP)=vecout(1)
!            RMSH(2,IP)=vecout(2)
!            RMSH(3,IP)=vecout(3)
!         ENDDO
!         print *,'yash123'
!         DO IC=MMSH1+1,MMSH1+NATMS
!            vecin(1)=RMSH(1,IC)
!            vecin(2)=RMSH(2,IC)
!            vecin(3)=RMSH(3,IC)
!            call rotvec(1,90,vecin,vecout)
!            RMSH(1,IC)=vecout(1)
!            RMSH(2,IC)=vecout(2)
!            RMSH(3,IC)=vecout(3)
!         ENDDO
! Checking non-orthonormality ERROR
!
         COULOMB(:) = 0.0D0
         RMSH = 0.0d0

         do i=1,NANG
         RMSH(1,i)=AMSH(1,i)
         RMSH(2,i)=AMSH(2,i)
         RMSH(3,i)=AMSH(3,i)
         enddo

         CALL COUPOT1

         TOLER = 1.0D-4
         ERROR = 1.0d-30

         DO I=1,NPOLY
         DO J=I,NPOLY
            SSSS=0.0D0
            DO IANG=1,NANG
               SSSS=SSSS+POLY(IANG,I)*POLY(IANG,J)*DOMEGANG(IANG)
            ENDDO
            IF(I.NE.J)THEN
            ERROR=MAX(ERROR,ABS(SSSS))
            IF((ABS(SSSS).GT.TOLER))PRINT *,I,J,SSSS
            ELSE
            !PRINT *,I,J,SSSS
            ENDIF
         ENDDO
         ENDDO

         SOLANG=0.0D0
         DO IANG=1,NANG
            SOLANG=SOLANG+DOMEGANG(IANG)
         ENDDO
         PRINT *,'SOLANG= ',SOLANG
         
         OPEN(199,FILE='SMULTI',form='formatted',status='unknown')
         REWIND(199)
!         icount2=0
         DO I=1,((LM+1)**2)
          SMULTI(I)=0.0D0
          kk=1
          DO IPTS=1,NANG
             ZPOTTEMP=0.0D0
             DO ICNT=1,NATMS
                DIST=0.0D0
                DO j=1,3
                DIST=DIST+(ATWD(j,ICNT)-AMSH(j,IPTS))**2
                ENDDO
                IF(DIST.gt.DELTA)THEN
                DIST=1.0D0/DSQRT(DIST)
                ZPOTTEMP=ZPOTTEMP-ZWATM(ICNT)*DIST
                ENDIF
             ENDDO
! COULOMB here is evaluated using the Wannier functions of central cell.            

             SMULTI(I)=SMULTI(I)+POLY(kk,I)*DOMEGANG(kk)
     &                *(COULOMB(IPTS)+ZPOTTEMP)
             kk=kk+1
          ENDDO ! IPTS
          ll=ceiling(dsqrt(dfloat(I)))
          SMULTI(I)=SMULTI(I)*(rad_sph**dfloat(ll))
!          IF(abs(SMULTI(I)).GT.1.0D-10)THEN
!          icount2=icount2+1
!          ENDIF
         ENDDO ! I

         WRITE(199,*)NPOLY  !icount2
           I=1
           SMULTI(I)=0.0d0
           WRITE(199,*)I,SMULTI(I) !zeroing monopole moment for a neutral system.
         DO I=2,((LM+1)**2)
!          IF(abs(SMULTI(I)).GT.1.0D-10)THEN
           WRITE(199,*)I,SMULTI(I)
!          ENDIF
         ENDDO
         CLOSE(199)

         CALL SYSTEM('cat SMULTI >> SMULTI_ALL')
         CALL SYSTEM('echo "===========================" >>SMULTI_ALL')
         CALL SYSTEM('echo "">> SMULTI_ALL')

! LATTICESUM: Evalautes contribution from outside Wannier Domain


!         CALL LATTICESUM
         ENDIF
!  777  continue

         FIRSTLAT=.TRUE.

         NGRAD=1
         IF ((IGGA(1).GT.0).OR.(IGGA(2).GT.0)) NGRAD=10

         do ipts=1,NMSH
            COULOMB(ipts)=COULCC(ipts)
            do IGRAD=1,NGRAD
               RHOG(ipts,IGRAD,1)=RHOCC(ipts,IGRAD,1)+RHOCC(ipts,1,2)
               RHOG(ipts,IGRAD,2)=RHOCC(ipts,IGRAD,2)
            enddo
         enddo


       ENDIF !FIRSTLAT

!===============================================================================
!===============================================================================

!      IF(IFOLD.EQ.2)THEN
!       PRINT *,'==============================================='
!       PRINT *,' UNFOLDING FROM CENTRAL CELL TO WANNIER DOMAIN'
!       PRINT *,'==============================================='
!
!       NMSW = NMSH
!       NMSH=NMSW*NDCELL
!
!
!       print *,'fold2',NMSW, NMSH
!
!       IF(NSPN.EQ.2)THEN
!       DO IPTS=1,NMSW
!          POT(IPTS)=POT(IPTS)
!          POT(IPTS+NMSH)=POT(IPTS+NMSW)
!       ENDDO
!       ENDIF
!
!       IFAK=NSPN-1
!       IOFS=IFAK*NMSH
!
!       NGRAD=1
!       IF ((IGGA(1).GT.0).OR.(IGGA(2).GT.0)) NGRAD=10
!       DO ILAT=1,NDCELL-1
!       DO IMSH=1,NMSW
!          POT(IMSH+ILAT*NMSW)=POT(IMSH)
!          POT(IMSH+IOFS+ILAT*NMSW)=POT(IMSH+IFAK*NMSH)
!          DO ISPN=1,NSPN
!          DO IGRAD=1,NGRAD
!             RHOG(IMSH+ILAT*NMSW,IGRAD,ISPN)=RHOG(IMSH,IGRAD,ISPN)
!          ENDDO
!          ENDDO
!       ENDDO
!       ENDDO
!
!       ENDIF

         DO i=1,NMSH
           RMSH(1,i)=RMSHCC(1,i)
           RMSH(2,i)=RMSHCC(2,i)
           RMSH(3,i)=RMSHCC(3,i)
         ENDDO
       RETURN
      END SUBROUTINE FOLDW
