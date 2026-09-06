c
c  routine to decompose the coulomb potential of the central cell into multipole moments (spherical harmonics)
c
       SUBROUTINE COUL_DECOMP
       INCLUDE 'PARAMS'
       INCLUDE 'commons.inc'
       PARAMETER (LM=06)
       PARAMETER (LMX=3*LM)
       PARAMETER (MAX_ANG=((LMX+1)*(LMX+2))/2)
       PARAMETER (MPX=MAX_ANG)
       LOGICAL ACTUAL
       COMMON/MIXPOT/POTIN(MAX_PTS*MXSPN),POTDV(MAX_PTS*MXSPN)
       COMMON/TMP1/POT(MAX_PTS),RHOG(MAX_PTS,10,MXSPN)
       DIMENSION ANG(3,MAX_ANG),DOMEGA(MAX_ANG)
       DIMENSION ZPOT(MAX_PTS*MXSPN)
       DIMENSION POLY(MPX,(LM+1)**2)
       DIMENSION SMULTI((LM+1)**2)
       DIMENSION ikk(49),CCENTER(3)
       DIMENSION atomvec(3),vec(3),amesh(3),atom1(3,100),zatm(100)

  
       ACTUAL = .TRUE.
       pi = 4.0D0*atan(1.0D0)

       print *, 'calling angmsh max_ang', MAX_ANG
c
c We get number  of angular grid points as NANG, grid point ANG, and DOMEGA as weight
c
        LMAX = LMX
        CALL ANGMSH(MAX_ANG,LMAX,NANG,ANG,DOMEGA)
        CALL flush(6)         
        NMSH = NANG
        print *, 'coul_decomp nmsh', nmsh
c
c choose radius and create grid on sphere. Increase the size of the 
c shpere and get new value of ANG as RMSH and WMSH as its corresponding weight.
c
         rad_sph = 30.0d0
         DO ipts = 1, NANG
           DO ix = 1,3
              rmsh(ix,ipts) = rad_sph*ang(ix,ipts)
           ENDDO
              wmsh(ipts) = domega(ipts)
         ENDDO

          do ipts=1,NANG
             print *,(rmsh(ix,ipts),ix=1,3)
          enddo       

         !CALL GCENTER(RMSH,NANG,CCENTER)


!         print *,'After shift'
!         do ipts=1,NANG
!            print *,(RMSH(ix,ipts),ix=1,3)
!         enddo       
!         stop

!        nmsh = nang
!        print *, 'nmsh, nang', nmsh, nang

! NANG is the total number of points on the sphere located at
! ANG(ix,NANG) locations such that HARMONICS gives Y_LM given by
! POLY(NANG,LM=J=1 to NPOLY (49)). For each grid point we have 49
! spherical harmonics.

         CALL HARMONICS(MPX,NANG,LM,ANG,POLY,NPOLY)
!        print *,NPOLY,(LM+1)**2
!        stop
!
! check for non-orthonormality ERROR
!
        TOLER = 1.0D-4
        ERROR = 1.0d-30

!        DO I=1,(LM+1)**2
!         DO J=I,(LM+1)**2
         DO I=1,NPOLY
         DO J=I,NPOLY
          SSSS=0.0D0
          DO IANG=1,NANG
           SSSS=SSSS+POLY(IANG,I)*POLY(IANG,J)*DOMEGA(IANG)
          END DO
          IF(I.NE.J)THEN
           ERROR=MAX(ERROR,ABS(SSSS))
           IF((ABS(SSSS).GT.TOLER))PRINT *,I,J,SSSS
          ELSE
           PRINT *,I,J,SSSS
          END IF
         END DO
         END DO
         SOLANG=0.0D0
         DO IANG=1,NANG
            SOLANG=SOLANG+DOMEGA(IANG)
         END DO
         PRINT *,'SOLANG= ',SOLANG

!       if(ACTUAL)then
!       CALL APOTNL(TOTQNUM,0,0)
!       do IANG=1,NMSH
!         x=RMSH(1,IANG)
!         y=RMSH(2,IANG)
!         z=RMSH(3,IANG)
!         r=dsqrt(x**2+y**2+z**2)
!!
!         POT(IANG) = (x**4+y**4+z**4 - 0.6d0*(x**2+y**2+z**2)**2)/r**4
!         !POT(IANG) = x*y*z/r**3 
!         !POT(IANG) = (x**6+y**6+z**6)/r**3 
!        enddo
       !CALL COUPOT1
!       COULOMB2=POT
!       else
!       CALL COULOMBPOT(RMSH,NANG,COULOMB2)
!        CALL COUPOT1
!       endif

!      print *,totqnum
!      print *,(POTDV(i),i=1,10)
!      print*,'------------'
!      print *,(POTIN(i),i=1,10)
!      print*,'------------'
!      print *,(COULOMB(i),i=1,10)



      open(99,file='XMOL.DAT',status='old')
      read(99,*)iatomtot
      read(99,*)

      atom1=0.0d0
      do i=1,iatomtot
       read(99,*)zatm(i),(atom1(ii,i),ii=1,3)
      enddo

      do igrid=1,NANG
        zpotemp=0.0d0
        amesh(:)=RMSH(:,igrid)
        do iatom=1,iatomtot
           atomvec(:)=atom1(:,iatom)
           vec(:)=atomvec(:)-amesh(:)
           denom=dsqrt(sum(vec(:)*vec(:)))
           zpotemp=zpotemp-zatm(iatom)/denom
        enddo
        ZPOT(igrid)=zpotemp 
      enddo

       OPEN(199,FILE='SMULTI',form='formatted',status='unknown')
       REWIND(199)
!       WRITE(199,*)rad_sph
!       OPEN(299,FILE='CPOTENTIAL',form='formatted',status='unknown')
!       REWIND(299)
!       OPEN(399,FILE='C2POTENTIAL',form='formatted',status='unknown')
!       REWIND(399)
!
!       DO IANG=1,NANG
!          WRITE(299,*)COULOMB2(IANG)
!       ENDDO

        icount2=0

        DO I=1,((LM+1)**2)
          SMULTI(I)=0.0D0
          DO IANG=1,NANG
!           SMULTI(I)=SMULTI(I)+POLY(IANG,I)*POTDV(IANG)*DOMEGA(IANG)
!           SMULTI(I)=SMULTI(I)+POLY(IANG,I)*
!     &        (POT(IANG)+ZPOT(IANG))*DOMEGA(IANG)
           SMULTI(I)=SMULTI(I)+POLY(IANG,I)*
     &        ZPOT(IANG)*DOMEGA(IANG)
!           SMULTI(I)=SMULTI(I)+POLY(IANG,I)*COULOMB2(IANG)*DOMEGA(IANG)
          END DO
          ll=ceiling(sqrt(dfloat(I))-1.0d0)
          fact1=dsqrt(4.0d0*pi/(2.0d0*ll+1.0d0))
          SMULTI(I)=SMULTI(I)*(rad_sph**dfloat(ll+1))/fact1
          smulti_toll=1.0D-15
          IF(abs(SMULTI(I)) .gt. smulti_toll)THEN
              icount2=icount2+1
          ENDIF
          WRITE(199,*)I,SMULTI(I)
        END DO

!        WRITE(199,*)(CCENTER(ii),ii=1,3)
!        WRITE(199,*)icount2
!        DO I=1,((LM+1)**2)
!         IF(abs(SMULTI(I)) .gt. smulti_toll)THEN
!         WRITE(199,*)I,SMULTI(I)
!         ENDIF 
!        ENDDO

!        DO IANG=1,NANG
!        potpot=0.0d0
!        DO I=1,49
!         ll=ceiling(sqrt(dfloat(I))-1)
!            potpot=potpot+SMULTI(I)*POLY(IANG,I)/rad_sph**(ll+1)
!        ENDDO
!        WRITE(399,*)potpot
!        ENDDO

! Stored SMULTI calculated considering rad_sph=5Bohr
!Reading SMULTI to calculate V(r) at a different sphere.
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!      
!       OPEN(150,FILE='SMULTI',form='formatted')
!       REWIND(150)
!
!       READ(150, *)rad_sph2
!       do I=1,48
!       READ(150, *)ikk(I),SMULTI(I)
!       enddo
!
!!       do ii=1,48
!!       print *,ikk(ii),SMULTI(ii)
!!       enddo
!
!       OPEN(151,FILE='testpot.dat',form='formatted')
!       REWIND(151)
!       WRITE(151,*)'old and new=',rad_sph2, rad_sph
!
!        DO IANG=1,NANG
!        potpot=0.0d0
!        DO I=1,49
!         ll=ceiling(sqrt(dfloat(I))-1)
!         potpot=potpot+(SMULTI(I)*POLY(IANG,I))/rad_sph**(ll+1)
!         ENDDO
!         WRITE(151,*)potpot
!        ENDDO
!
!        close(150)
!        close(151)
!!       stop
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

        CLOSE(199)
!        CLOSE(299)
!        CLOSE(399)
        call stopit

      end subroutine coul_decomp

         
  



