      SUBROUTINE FOLDW(IFOLD)
      INCLUDE 'PARAMS'
      INCLUDE 'commons.inc'
      PARAMETER (LM=06)
      PARAMETER (LMX=3*LM)
      PARAMETER (MAX_ANG=((LMX+1)*(LMX+2))/2)
      PARAMETER (MPX=MAX_ANG)
      PARAMETER (MXPR=MXPOISS)
      DIMENSION POLY(MPX,(LM+1)**2),SMULTI((LM+1)**2)
      DIMENSION vecin(3),vecout(3)
      LOGICAL FIRST,FIRSTLAT
      COMMON/FLOMESH/RABCD(3,1000)
      COMMON/MIXPOT/POTIN(MAX_PTS*MXSPN),POT(MAX_PTS*MXSPN)
      COMMON/TMP1/COULOMB(MAX_PTS),RHOG(MAX_PTS,NVGRAD,MXSPN)
      COMMON/REPLICA/NMSHW,NDCELL,MMSH1,NATMS,
     &      ZWATM(1000),RATMS(3,1000),ITSCF
      COMMON/PTRANS/TVEC(3,3),ATHETA,IRAXIS,NXYZ,iwx,iwy,iwz
      COMMON/COULDW/NANG,ANG(3,MAX_ANG),DOMEGANG(MAX_ANG),RAD_SPH
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
!        open(402,file="RHOWC")
!        rewind(402)
!        DO IMSH=1,NMSHW
!        DO ISPN=1,NSPN
!        DO IGRAD=1,NGRAD
!!           WRITE(402,*)RHOG(IMSH,IGRAD,ISPN),WMSH(IMSH)
!          READ(402,*)RHOG(IMSH,IGRAD,ISPN)!,WMSH(IMSH)
!        ENDDO
!        ENDDO
!        ENDDO 
!        close(402)
!================================================================================
! V_ext(eN) from Wannier Domain evaluated in Central Cell.
       DELTA=1.0D-100
       POTIN=0.0D0
       DO iatm=MMSH1+1,MMSH1+NDCELL*NATMS
          vx=RMSH(1,iatm)
          vy=RMSH(2,iatm)
          vz=RMSH(3,iatm)
       DO IMSH=1,MMSH1
          DIST = (RMSH(1,IMSH)-vx)**2+(RMSH(2,IMSH)-vy)**2+
     &           (RMSH(3,IMSH)-vz)**2
          DIST=MAX(DIST,DELTA)
          DIST=1.0D0/DSQRT(DIST)
          POTIN(IMSH)= POTIN(IMSH)-ZWATM(iatm-MMSH1)*DIST !make ZATMS general
       ENDDO
       ENDDO

       CALL GTTIME(TPOT2)

       PRINT *,'TIME ELAPSED IN POTINW',(TPOT2-TPOT1)

       PRINT *,'==============================================='
       PRINT *,'  FOLDING FROM WANNIER DOMAIN TO CENTRAL CELL'
       PRINT *,'==============================================='

! CENTRAL CELL CHARGE BEFORE FOLDING

        qtest2=0.0d0
        qtest=0.0d0
        DO ipts=1,NMSHW
          qtest=qtest+RHOG(ipts,1,1)*WMSH(ipts)
        ENDDO
        PRINT *,'central_cell',qtest

        qtest2=qtest
        ilatcell=-(NDCELL-1)/2-1
        DO ILAT=1,NDCELL-1
        qtest1=0.0d0
        ilatcell=ilatcell+1
        if(ilatcell.eq.0)ilatcell=ilatcell+1
        DO ipts=1,NMSHW
          qtest1=qtest1+RHOG(ipts+ILAT*NMSHW,1,1)*WMSH(ipts+ILAT*NMSHW)
        ENDDO
        qtest2=qtest2+qtest1
        print *,'CELL qtest1',ilatcell,qtest1
        ENDDO

        print *,'Total Charge:',qtest2

! FOLDING COULOMB POTENTIAL TO CENTRAL CELL

        DO ILAT=1,NDCELL-1
        DO IMSH=1,NMSHW
           COULOMB(IMSH)=COULOMB(IMSH)+COULOMB(IMSH+ILAT*NMSHW)
        ENDDO 
        ENDDO

        DO IMSH=1,NMSHW
           COULOMB(IMSH)=COULOMB(IMSH)+POTIN(IMSH)
        ENDDO 

! FOLDING WANNIER CHARGE DENSITY TO CENTRAL CELL
        NGRAD=1
        IF ((IGGA(1).GT.0).OR.(IGGA(2).GT.0)) NGRAD=10
        DO ILAT=1,NDCELL-1
        DO IMSH=1,NMSHW
           DO ISPN=1,NSPN
           DO IGRAD=1,NGRAD
              RHOG(IMSH,IGRAD,ISPN)=RHOG(IMSH,IGRAD,ISPN)+
     &            RHOG(IMSH+ILAT*NMSHW,IGRAD,ISPN)
           ENDDO
           ENDDO
        ENDDO 
        ENDDO

! TOTAL CHARGE AFTER FOLDING 
        qtest2=0.0d0
        qtest3=0.0d0
        DO ipts=1,NMSHW
          qtest2=qtest2+RHOG(ipts,1,1)*WMSH(ipts)
          qtest3=qtest3+abs((RHOG(ipts,1,1)-2.0d0*RHOG(ipts,1,2)))
     &           *WMSH(ipts)
        ENDDO

       PRINT *,'delta_rho, CC(q), FOLD=',qtest3,qtest,qtest2

        ECOUL=0.0D0
        DO IPTS=1,NMSHW
           ECOUL = ECOUL+0.5D0*COULOMB(IPTS)*RHOG(IPTS,1,1)*WMSH(IPTS)
        ENDDO 

!==========================================================================
       PRINT *,'NATMS AND Z VALUE',NATMS,(ZWATM(ii),ii=1,NATMS)
       DO ii=MMSH1+1,MMSH1+NATMS
          PRINT *,(RMSH(jj,ii),jj=1,3)
       ENDDO

!==========================================================================
! Coulomb potential is evaluated at each nuclear site in the Wannier domain 
! and then integrated
!==========================================================================
       POTIN=0.0d0
       DO iat=1,NATMS
       DO iwd=0,NDCELL-1
          pottemp=COULOMB(MMSH1+iat+iwd*natms)
          POTIN(iat)=POTIN(iat)+pottemp
       ENDDO
       ENDDO

       ELOCALW=0.0D0
       DO j=1,NATMS
       ELOCALW=ELOCALW-0.5D0*ZWATM(j)*POTIN(j)
       ENDDO
!==========================================================================
! ENNUCW
!==========================================================================
       IF(FIRST)THEN
         PRINT *,'Nuclear coordinates in Wannier Domain'
         DO ii=MMSH1+1,MMSH1+NATMS*NDCELL
         PRINT *,(RMSH(j,ii),j=1,3),ZWATM(ii-MMSH1)
         ENDDO

         POTIN=0.0d0
         DO iat=1,NATMS
         DO jat=1,NATMS
         DO iwd=0,NDCELL-1
         vx=RMSH(1,MMSH1+iwd*natms+jat)
         vy=RMSH(2,MMSH1+iwd*natms+jat)
         vz=RMSH(3,MMSH1+iwd*natms+jat)

         DIST = (RMSH(1,MMSH1+iat)-vx)**2+
     &          (RMSH(2,MMSH1+iat)-vy)**2+
     &          (RMSH(3,MMSH1+iat)-vz)**2
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
          DO IPTS=MMSH1+NDCELL*NATMS+1,NMSH
             ZPOTTEMP=0.0D0
             DO ICNT=MMSH1+1,MMSH1+NATMS
                DIST=0.0D0
                DO j=1,3
                DIST=DIST+(RMSH(j,ICNT)-RMSH(j,IPTS))**2
                ENDDO
                DIST=1.0D0/DSQRT(DIST)
                ZPOTTEMP=ZPOTTEMP-ZWATM(ICNT-MMSH1)*DIST
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
         !CALL LATTICESUM
         ENDIF
!  777  continue

         FIRSTLAT=.TRUE.
         NMSH=NMSHW
       ENDIF !FIRSTLAT

!===============================================================================
!===============================================================================

      IF(IFOLD.EQ.2)THEN
       PRINT *,'==============================================='
       PRINT *,' UNFOLDING FROM CENTRAL CELL TO WANNIER DOMAIN'
       PRINT *,'==============================================='

       NMSH=MMSH1+NDCELL*NATMS+NANG


       IF(NSPN.EQ.2)THEN
       DO IPTS=1,NMSHW
          POT(IPTS)=POT(IPTS)
          POT(IPTS+NMSH)=POT(IPTS+NMSHW)
       ENDDO
       ENDIF
      
       IFAK=NSPN-1
       IOFS=IFAK*NMSH

       NGRAD=1
       IF ((IGGA(1).GT.0).OR.(IGGA(2).GT.0)) NGRAD=10
       DO ILAT=1,NDCELL-1
       DO IMSH=1,NMSHW
          POT(IMSH+ILAT*NMSHW)=POT(IMSH)
          POT(IMSH+IOFS+ILAT*NMSHW)=POT(IMSH+IFAK*NMSH)
         ! DO ISPN=1,NSPN
         ! DO IGRAD=1,NGRAD
         !    RHOG(IMSH+ILAT*NMSHW,IGRAD,ISPN)=RHOG(IMSH,IGRAD,ISPN)
         ! ENDDO
         ! ENDDO
       ENDDO 
       ENDDO

       ENDIF

       RETURN
      END SUBROUTINE FOLDW
