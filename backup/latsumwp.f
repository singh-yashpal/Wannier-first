      SUBROUTINE LATTICESUM
      INCLUDE 'PARAMS'
      INCLUDE 'commons.inc'
      INCLUDE 'mpif.h'
      DIMENSION rr1(3),coord(3),rcellw(3,2000),icellw(3,2000)
      DIMENSION atoms(3,1000)
      LOGICAL FGRID
      COMMON/MIXPOT/POTIN(MAX_PTS*MXSPN),POT(MAX_PTS*MXSPN)
      COMMON/TMP1/COULOMB(MAX_PTS),RHOG(MAX_PTS,NVGRAD,MXSPN)
      COMMON/REPLICA/NMSHW,NDCELL,MMSH1,NATMS,
     &       ZWATM(1000),RATMS(3,1000),ITSCF
      COMMON/PTRANS/TVEC(3,3),atheta,iraxis,NXYZ,iwx,iwy,iwz
      COMMON/LMULTI/ILMTOTW,ILMW(100),SMULTIW(100)
      DATA FGRID/.TRUE./

      NMSHW=MMSH1/NDCELL
      pi = 4.0D0*atan(1.0D0)

      PRINT *,'NMSHW, MMSH1, NDCELL',NMSHW,MMSH1,NDCELL
! To remove the lattice points already considered in the calculation.
      open(20,file='UTRANS')
      rewind(20)
      rcellw=0.0d0
      do i=1,1000
      read(20,*,end=1011)rcellw(1,i),rcellw(2,i),rcellw(3,i)
      enddo
 1011 continue
      close(20)


      OPEN(299,FILE='SMULTI',FORM='formatted',STATUS='OLD')
      READ(299,*)ILMTOTW
      DO KK=1,ILMTOTW
         READ(299,*)ILMW(KK),SMULTIW(KK)
      ENDDO
      CLOSE(299)
!============================================================
      PRINT *,'RUNNING LATTICE SUM, NXYZ=',NXYZ
      IF(FGRID)THEN
      ncellw=NDCELL-1
      do kk=1,ncellw
      print *,'cells',(rcellw(j,kk),j=1,3)
      enddo

      do ii=1,3
      print *,'TVEC',(TVEC(jj,ii),jj=1,3)
      enddo
      PRINT *,'CELLS IN WANNIER DOMAIN',ncellw
!============================================================
      IF(ncellw.GT.0)THEN 
         PRINT *,'FGRID IN LATTICESUM'
         PRINT *,'Calling LSGRID'
!         CALL LSGRID(rcellw,ncellw,nxyz)
         CALL LSGRIDD(rcellw,ncellw) !Evaluating grid points for the direct calculation 
!         CALL LSGRIDW(rcellw,ncellw) !Grid from gaussp
         FGRID=.FALSE.
         PRINT *,'LSGRID DONE'
      ENDIF
      ENDIF
!============================================================

! LATSD: Purely for 1 dimensional direct calculation.
! LAT_DIRECT: Direct calculation for 1,2 or 3D       

      CALL GTTIME(TLAT1)
      ELATSUM = 0.0D0

      IF (NPROC.GT.0) THEN
         POT = 0.0D0
         CALL SENDDATA(108)
         CALL PAMLATSUM(NMSHW)
         DO IPTS=1,NMSHW
            ELATSUM = ELATSUM+0.5D0*POT(IPTS)
     &               *RHOG(IPTS,1,1)*WMSH(IPTS)
            COULOMB(IPTS)=COULOMB(IPTS)+POT(IPTS)
         ENDDO
      ELSE
         DO IPTS = 1, NMSHW
            PHIW = 0.0D0
            DO IX=1,3
               RR1(IX)=RMSH(IX,IPTS)
            ENDDO
            CALL LAT_DIRECT(RR1,PHIW)
            ELATSUM = ELATSUM+0.5D0*PHIW*RHOG(IPTS,1,1)*WMSH(IPTS)
            COULOMB(IPTS)=COULOMB(IPTS)+PHIW
         ENDDO
      ENDIF


      EZNUCLAT = 0.0D0
      DO IPTS = 1, NATMS
         DO IX=1,3
            RR1(IX)=RMSH(IX,MMSH1+IPTS)
         ENDDO
         PHIW = 0.0D0
         CALL LAT_DIRECT(RR1,PHIW)
         EZNUCLAT=EZNUCLAT-0.5D0*PHIW*ZWATM(IPTS)
      ENDDO





!       COULOMBW=0.0D0
!       DO ipts=1,NMSHW
!          phiw=0.0d0
!          DO ix=1,3
!             rr1(ix)=RMSH(ix,ipts)
!          ENDDO
!          CALL LAT_DIRECT(rr1,phiw)
!!          IF(NXYZ.eq.1.and.ncellw.gt.0)then
!!          !CALL LATSD(TVEC,rr1,rcellw,ncellw,phiw)
!!          CALL LAT_SUM(rr1,phiw)
!!          ENDIF
!!
!!          IF(NXYZ.eq.2.and.ncellw.gt.0)then
!!          CALL LAT_SUM2(TVEC,rr0,atoms,natoms,phiw)
!!          ENDIF
!!          IF(NXYZ.eq.3.and.ncellw.gt.0)then
!!          CALL LAT_SUM(TVEC,rr0,rcellw,ncellw,phiw)
!!          ENDIF
!          COULOMBW(ipts)=phiw
!       ENDDO
!
!       EZNUCLAT=0.0D0
!       DO ii=1,NATMS
!          DO ix=1,3
!             rr1(ix)=RMSH(ix,MMSH1+ii)
!          ENDDO
!          phiw=0.0d0
!          CALL LAT_DIRECT(rr1,phiw)
!!          IF(NXYZ.eq.1 .and. ncellw.gt.0)then
!!          CALL LATSD(TVEC,rr1,rcellw,ncellw,phiw)
!!          CALL LAT_SUM1(rr1,phiw)
!!          ENDIF
!
!!          IF(NXYZ.eq.2.and.ncellw.gt.0)then
!!          CALL LAT_SUM2(TVEC,rr0,atoms,natoms,phiw)
!!          CALL LAT_DIRECT(rr1,phiw)
!!          ENDIF
!
!!          IF(NXYZ.eq.3.and.ncellw.gt.0)then
!!          CALL LAT_SUM(TVEC,rr0,rcellw,ncellw,phiw)
!!          CALL LAT_DIRECT(rr1,phiw)
!!          ENDIF
!          EZNUCLAT=EZNUCLAT-0.5D0*phiw*ZWATM(ii)
!       ENDDO

      CALL GTTIME(TLAT2)
      PRINT *,'TLATSUM',(TLAT2-TLAT1)
       PRINT *,'EZNUCLAT=',EZNUCLAT
       PRINT *,'END LATTICE SUM'
       RETURN
      END SUBROUTINE LATTICESUM

!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      SUBROUTINE LSGRIDD(rcellw,ncellw)
      IMPLICIT REAL*8 (a-h,o-z)
      PARAMETER (ILAT=20)
      PARAMETER (IGRID=(2*ILAT+1)**3)
      DIMENSION rcellw(3,ncellw),ngrid(3,IGRID),angrid(IGRID)
      COMMON/PTRANS/TVEC(3,3),atheta,iraxis,
     &             NXYZ,ISHELLX,ISHELLY,ISHELLZ

      eps=1.0d-2
      iwx=0
      iwy=0
      iwz=0

      if (NXYZ.ge.1)iwx=ILAT
      if (NXYZ.ge.2)iwy=ILAT
      if (NXYZ.ge.3)iwz=ILAT

C From negative to positive values of n1,n2, and n3
      icount=0
      do n1=-iwx,iwx
      do n2=-iwy,iwy
      do n3=-iwz,iwz

!=================================================================
         if(n1.eq.0.and.n2.eq.0.and.n3.eq.0)goto 800
! Removing cells considered in Wannier Domain
!=================================================================
         xn1=dfloat(n1)
         xn2=dfloat(n2)
         xn3=dfloat(n3)

         do jj=1,ncellw
         if(abs(xn1-rcellw(1,jj)).lt.eps .and. abs(xn2-rcellw(2,jj))
     &     .lt.eps .and. abs(xn3-rcellw(3,jj)).lt.eps)goto 800
         enddo
!=================================================================
         icount=icount+1
         ngrid(1,icount)=n1
         ngrid(2,icount)=n2
         ngrid(3,icount)=n3
         if(iraxis.eq.1)then
         angrid(icount)=atheta*dfloat(n1)
         endif
         if(iraxis.eq.2)then
         angrid(icount)=atheta*dfloat(n2)
         endif
         if(iraxis.eq.3)then
         angrid(icount)=atheta*dfloat(n3)
         endif
  800 continue
      enddo
      enddo
      enddo

      open(30,file="NGRID",form="formatted")
      rewind(30)
      WRITE(30,*)icount
      DO ii=1,icount
      WRITE(30,'(3(1X,I3),(1X,f12.6))')(ngrid(ix,ii),ix=1,3),angrid(ii)
      ENDDO

      print *,'Total gird points outside WD:',icount
      close(30)

      RETURN
      END SUBROUTINE LSGRIDD

!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      SUBROUTINE LSGRIDW(rcellw,ncellw,nxyz)
      IMPLICIT REAL*8 (a-h,o-z)
      PARAMETER (nn=5) !Gauss quadrature nn points
      PARAMETER (l=10) !Direct sum part

      PARAMETER (igrid=(l**3+3*nn*(l**2)+6*(nn**2)*l+6*nn**3))
      PARAMETER (iygrid=(8*igrid))
      DIMENSION sump(10),alp(10)
      DIMENSION xgrid(3,igrid),wgrid(3,igrid)
      DIMENSION ygrid(3,iygrid),wygrid(3,iygrid)

      DIMENSION xn1(10),xn2(10),xn3(10),wn1(10),wn2(10),wn3(10)
!      COMMON/LATGRID/XVEC(10,MAX_PTS),YVEC(10,MAX_PTS),ZVEC(10,MAX_PTS),
!     &    WXVEC(10,MAX_PTS),WYVEC(10,MAX_PTS),WZVEC(10,MAX_PTS),IGRID
      DIMENSION rcellw(3,ncellw),ngrid(3,IGRID),angrid(IGRID)

      Character(8) :: Date
      Character(10) :: Time1,Time2
      Character(5) :: Zone
      Integer :: Values(8)

      pi=4.0D0*atan(1.0D0)

      print *,'nn,l,nxyz',nn,l,nxyz

      do m=1,10
         sump(m)=0.0D0
         alp(m)=0.0D0
      enddo

      DO n=1,l-1
      DO m=2,10
      sump(m)=sump(m)+1.0D0/(dfloat(n)**dfloat(m))
      ENDDO
      ENDDO

      alp(2)=pi*pi/6.0D0-sump(2)
      alp(2)=1.0D0/alp(2)
      alp(3)=2.0D0*(1.20205690315959D0-sump(3))
      alp(3)=1.0D0/alp(3)**(1.0D0/2.0D0)
      alp(4)=3.0D0*(1.08232323371114D0-sump(4))
      alp(4)=1.0D0/alp(4)**(1.0D0/3.0D0)
      alp(5)=alp(4)
! Zeta(5) is blowing and for now we approximate zeta(5) with zeta(4)
!      alp(5)=4.0D0*(1.0369277551433699263-sump(5))
!      alp(5)=1.0D0/alp(5)**(1.0D0/4.0D0)
      alp(6)=5.0D0*(pi**6.0d0/945.0D0-sump(6))
      alp(6)=1.0D0/alp(6)**(1.0D0/5.0D0)

      DO ii=7,10
      alp(ii)=alp(6)
      ENDDO

      xgrid=0.0d0
      wgrid=0.0d0

!      OPEN(299,FILE='SMULTI',FORM='formatted',STATUS='OLD')
!      REWIND(299)
!      READ(299,*)ILMTOT
!      m_old=0
!
!      DO 505 kk=1,ILMTOT
!      READ(299,*)ILM(kk),SMULTI(kk)
!
!      mm=ILM(kk) 
!      m=ceiling(dsqrt(dfloat(mm)))
!      if(m.LE.m_old)goto 505
!
      icount=0
      m=2

      iwx=1
      iwy=1
      iwz=1

      if(NXYZ.ge.1)iwx=l
      if(NXYZ.ge.2)iwy=l
      if(NXYZ.ge.3)iwz=l

      do n3=0,iwz-1
      do n2=0,iwy-1
      do n1=0,iwx-1

         print *,'gridd',n1,n2,n3
         if(n1.eq.0.and.n2.eq.0.and.n3.eq.0)goto 800
! Removing cells considered in Wannier Domain
!=================================================================
         do jj=1,ncellw
         if(n1.eq.rcellw(1,jj).and.n2.eq.rcellw(2,jj).and.
     &      n3.eq.rcellw(3,jj))goto 800
         enddo
!=================================================================
         icount=icount+1
         xgrid(1,icount)=dfloat(n1)
         xgrid(2,icount)=dfloat(n2)
         xgrid(3,icount)=dfloat(n3)

         wgrid(1,icount)=1.0d0
         wgrid(2,icount)=1.0d0
         wgrid(3,icount)=1.0d0

  800    continue
         enddo ! n1 loop ends
!==========================================================================
         dlimn10 = 0.0D0
         dlimn11 = 1.0D0/alp(m)

         call gaussp(dlimn10,dlimn11,nn,xn1,wn1)
         do i=1,nn
            icount=icount+1
            xgrid(1,icount)= 1.0d0/xn1(i)
            xgrid(2,icount)= dfloat(n2)
            xgrid(3,icount)= dfloat(n3)

            wgrid(1,icount)= wn1(i)/(xn1(i)**2.0d0)
            wgrid(2,icount)= 1.0d0
            wgrid(3,icount)= 1.0d0
         enddo

         if(nxyz.eq.1)goto 601
       enddo !n2 loop ends      
!------------------------------------------------------------------------            
       dlimn20=0.0D0
       dlimn21=pi/4.0D0
       call gaussp(dlimn20,dlimn21,nn,xn2,wn2)
       do j=1,nn
          dlimn10=0.0D0
          dlimn11=1.0D0/(alp(m)*cos(xn2(j)))
          call gaussp(dlimn10,dlimn11,nn,xn1,wn1)
          do i=1,nn

            rr=xn1(i)
            vv=rr*cos(xn2(j))
            ww=rr*sin(xn2(j))

            icount=icount+1
            xgrid(3,icount)=dfloat(n3)
            xgrid(2,icount)=1.0d0/vv
            xgrid(1,icount)=1.0d0/ww

            wgrid(3,icount)=1.0d0
            wgrid(2,icount)=wn2(j)*(dsqrt(rr)/vv**2.0d0)
            wgrid(1,icount)=wn1(i)*(dsqrt(rr)/ww**2.0d0)

!======================================================================              
!           Term 2
!======================================================================              
            vv=rr*sin(xn2(j))
            ww=rr*cos(xn2(j))

            icount=icount+1
            xgrid(3,icount)=dfloat(n3)
            xgrid(2,icount)=1.0d0/vv
            xgrid(1,icount)=1.0d0/ww

            wgrid(3,icount)=1.0d0
            wgrid(2,icount)=wn2(j)*(dsqrt(rr)/vv**2.0d0)
            wgrid(1,icount)=wn1(i)*(dsqrt(rr)/ww**2.0d0)
         enddo
      enddo
 503  continue
      enddo ! n3 loop ends.       
!=====================================================================              
!      do n3=0,l-1
!        do n1=0,l-1
!          dlimn20 = 0.0D0
!          dlimn21 = 1.0D0/alp(m)
!          call gaussp(dlimn20,dlimn21,nn,xn2,wn2)
!          do j=1,nn
!            iicount=iicount+1
!
!            wxvec(m,iicount)=1.0d0
!            wyvec(m,iicount)=wn2(j)/(xn2(j)**2.0d0)
!            wzvec(m,iicount)=1.0d0
!            xvec(m,iicount)=dfloat(n1)
!            yvec(m,iicount)=1.0d0/xn2(j)
!            zvec(m,iicount)=dfloat(n3)
!          enddo
!         enddo ! n1 loop ends
!!=======================================================================
!
!         dlimn20 = 0.0D0
!         dlimn21 = pi/4.0D0
!         call gaussp(dlimn20,dlimn21,nn,xn2,wn2)
!         do j=1,nn
!            dlimn10 = 0.0D0
!            dlimn11 = 1.0D0/alp(m)*cos(xn2(j))
!            call gaussp(dlimn10,dlimn11,nn,xn1,wn1)
!            do i=1,nn
!              rr=xn1(i)
!              uu=rr*cos(xn2(j))
!              vv=rr*sin(xn2(j))
!              iicount=iicount+1
!
!              wxvec(m,iicount)=wn1(i)*(dsqrt(rr)/uu**2.0d0)
!              wyvec(m,iicount)=wn2(j)*(dsqrt(rr)/vv**2.0d0)
!              wzvec(m,iicount)=1.0d0
!              xvec(m,iicount)=1.0d0/uu
!              yvec(m,iicount)=1.0d0/vv
!              zvec(m,iicount)=dfloat(n3)
!
!!=====================================================================
!              ! PART 2
!
!              uu=rr*sin(xn2(j))
!              vv=rr*cos(xn2(j))
!
!              iicount=iicount+1
!
!              wxvec(m,iicount)=wn1(i)*(dsqrt(rr)/uu**2.0d0)
!              wyvec(m,iicount)=wn2(j)*(dsqrt(rr)/vv**2.0d0)
!              wzvec(m,iicount)=1.0d0
!              xvec(m,iicount)=1.0d0/uu
!              yvec(m,iicount)=1.0d0/vv
!              zvec(m,iicount)=dfloat(n3)
!
!            enddo
!         enddo
!        enddo ! n3 loop      
!!========================================================================
!
!        do n2=0,l-1
!           do n3=0,l-1
!              dlimn10 = 0.0D0
!              dlimn11 = 1.0D0/alp(m)
!              call gaussp(dlimn10,dlimn11,nn,xn1,wn1)
!              do i=1,nn
!                 iicount=iicount+1
!
!                 wxvec(m,iicount)=wn1(i)/(xn1(i)**2.0d0)
!                 wyvec(m,iicount)=1.0d0
!                 wzvec(m,iicount)=1.0d0
!                 xvec(m,iicount)=1.0d0/xn1(i)
!                 yvec(m,iicount)=dfloat(n2)
!                 zvec(m,iicount)=dfloat(n3)
!               enddo
!           enddo ! n3 loop ends
!
!!====================================================================================
!
!           dlimn10 = 0.0D0
!           dlimn11 = pi/4.0D0
!           call gaussp(dlimn10,dlimn11,nn,xn1,wn1)
!           do i=1,nn
!              dlimn30 = 0.0D0
!              dlimn31 = 1.0D0/(alp(m)*cos(xn1(i)))
!              call gaussp(dlimn30,dlimn31,nn,xn3,wn3)
!              do k=1,nn
!                 rr=xn3(k)
!                 uu=rr*cos(xn1(i))
!                 ww=rr*sin(xn1(i))
!
!                 iicount=iicount+1
!
!                 wxvec(m,iicount)=wn1(i)*(dsqrt(rr)/uu**2.0d0)
!                 wyvec(m,iicount)=1.0d0
!                 wzvec(m,iicount)=wn3(k)*(dsqrt(rr)/ww**2.0d0)
!                 xvec(m,iicount)=1.0d0/uu
!                 yvec(m,iicount)=dfloat(n2)
!                 zvec(m,iicount)=1.0d0/ww
!
!!=========================================================================
!                 uu=rr*sin(xn1(i))
!                 ww=rr*cos(xn1(i))
!
!                 iicount=iicount+1
!
!                 wxvec(m,iicount)=wn1(i)*(dsqrt(rr)/uu**2.0d0)
!                 wyvec(m,iicount)=1.0d0
!                 wzvec(m,iicount)=wn3(k)*(dsqrt(rr)/ww**2.0d0)
!                 xvec(m,iicount)=1.0d0/uu
!                 yvec(m,iicount)=dfloat(n2)
!                 zvec(m,iicount)=1.0d0/ww
!               enddo
!            enddo
!        enddo !n2 loop
!
!!======================================================================
!        ! Triple Integration
!!======================================================================
!        dlimn30 = 0.0D0
!        dlimn31 = pi/4.0D0
!        call gaussp(dlimn30,dlimn31,nn,xn3,wn3)
!        do k=1,nn
!           dlimn20 = atan(1.0D0/sin(xn3(k)))
!           dlimn21 = pi/2.0D0
!           call gaussp(dlimn20,dlimn21,nn,xn2,wn2)
!           do j=1,nn
!              dlimn10 = 0.0D0
!              dlimn11 = 1.0D0/(alp(m)*sin(xn2(j))*cos(xn3(k)))
!              call gaussp(dlimn10,dlimn11,nn,xn1,wn1)
!              do i=1,nn
!! Term 1
!                 rr=xn1(i)
!                 wttfact=(rr*rr*sin(xn2(j)))**(1.0d0/3.0d0)
!  
!                 uu=rr*sin(xn2(j))*cos(xn3(k))
!                 vv=rr*sin(xn2(j))*sin(xn3(k))
!                 ww=rr*cos(xn2(j))
!  
!                 iicount=iicount+1
!
!                 wxvec(m,iicount)=(wn1(i)*wttfact)/(uu**2.0d0)
!                 wyvec(m,iicount)=(wn2(j)*wttfact)/(vv**2.0d0)
!                 wzvec(m,iicount)=(wn3(k)*wttfact)/(ww**2.0d0)
!                 xvec(m,iicount)=1.0d0/uu
!                 yvec(m,iicount)=1.0d0/vv
!                 zvec(m,iicount)=1.0d0/ww
!! Term 2
!                 uu=rr*sin(xn2(j))*sin(xn3(k))
!                 vv=rr*sin(xn2(j))*cos(xn3(k))
!                 ww=rr*cos(xn2(j))
!  
!                 iicount=iicount+1
!
!                 wxvec(m,iicount)=(wn1(i)*wttfact)/(uu**2.0d0)
!                 wyvec(m,iicount)=(wn2(j)*wttfact)/(vv**2.0d0)
!                 wzvec(m,iicount)=(wn3(k)*wttfact)/(ww**2.0d0)
!                 xvec(m,iicount)=1.0d0/uu
!                 yvec(m,iicount)=1.0d0/vv
!                 zvec(m,iicount)=1.0d0/ww
!! Term   3
!                 uu=rr*sin(xn2(j))*sin(xn3(k))
!                 vv=rr*cos(xn2(j))
!                 ww=rr*sin(xn2(j))*cos(xn3(k))
!   
!                 iicount=iicount+1
!
!                 wxvec(m,iicount)=(wn1(i)*wttfact)/(uu**2.0d0)
!                 wyvec(m,iicount)=(wn2(j)*wttfact)/(vv**2.0d0)
!                 wzvec(m,iicount)=(wn3(k)*wttfact)/(ww**2.0d0)
!                 xvec(m,iicount)=1.0d0/uu
!                 yvec(m,iicount)=1.0d0/vv
!                 zvec(m,iicount)=1.0d0/ww
!! Term 4
!                 uu=rr*sin(xn2(j))*cos(xn3(k))
!                 vv=rr*cos(xn2(j))
!                 ww=rr*sin(xn2(j))*sin(xn3(k))
!
!                 iicount=iicount+1
!
!                 wxvec(m,iicount)=(wn1(i)*wttfact)/(uu**2.0d0)
!                 wyvec(m,iicount)=(wn2(j)*wttfact)/(vv**2.0d0)
!                 wzvec(m,iicount)=(wn3(k)*wttfact)/(ww**2.0d0)
!                 xvec(m,iicount)=1.0d0/uu
!                 yvec(m,iicount)=1.0d0/vv
!                 zvec(m,iicount)=1.0d0/ww
!! Term 5
!                 xx=rr*cos(xn2(j))
!                 vv=rr*sin(xn2(j))*cos(xn3(k))
!                 ww=rr*sin(xn2(j))*sin(xn3(k))
!
!                 iicount=iicount+1
!
!                 wxvec(m,iicount)=(wn1(i)*wttfact)/(uu**2.0d0)
!                 wyvec(m,iicount)=(wn2(j)*wttfact)/(vv**2.0d0)
!                 wzvec(m,iicount)=(wn3(k)*wttfact)/(ww**2.0d0)
!                 xvec(m,iicount)=1.0d0/uu
!                 yvec(m,iicount)=1.0d0/vv
!                 zvec(m,iicount)=1.0d0/ww
!! Term 6
!                 xx=rr*cos(xn2(j))
!                 vv=rr*sin(xn2(j))*sin(xn3(k))
!                 ww=rr*sin(xn2(j))*cos(xn3(k))
!   
!                 iicount=iicount+1
!
!                 wxvec(m,iicount)=(wn1(i)*wttfact)/(uu**2.0d0)
!                 wyvec(m,iicount)=(wn2(j)*wttfact)/(vv**2.0d0)
!                 wzvec(m,iicount)=(wn3(k)*wttfact)/(ww**2.0d0)
!                 xvec(m,iicount)=1.0d0/uu
!                 yvec(m,iicount)=1.0d0/vv
!                 zvec(m,iicount)=1.0d0/ww
!          enddo
!        enddo       
!       enddo       
!
!      m_old=m

 601    CONTINUE ! m loop started at the top ends here 

!      CLOSE(299)

! Copying grid points in all the quadrants.
!  1  1  1
!  1  1 -1
!  1 -1  1
! -1  1  1
!  1 -1 -1
! -1  1 -1
! -1 -1  1
! -1 -1 -1


      jcount=0

      do ii=1,icount

      xn=xgrid(1,ii)
      yn=xgrid(2,ii)
      zn=xgrid(3,ii)

      wx=wgrid(1,ii)
      wy=wgrid(2,ii)
      wz=wgrid(3,ii)
      
!1------------- 1 0 0------------------------------
          if(xn.ne.0.0d0.and.yn.eq.0.0d0.and.zn.eq.0.0d0)then
          jcount=jcount+1
          ygrid(1,jcount)=xn
          ygrid(2,jcount)=yn
          ygrid(3,jcount)=zn

          wygrid(1,jcount)=wx
          wygrid(2,jcount)=wy
          wygrid(3,jcount)=wz


!1------------- -1 0 0------------------------------
          jcount=jcount+1
          ygrid(1,jcount)=-xn
          ygrid(2,jcount)=yn
          ygrid(3,jcount)=zn

          wygrid(1,jcount)=wx
          wygrid(2,jcount)=wy
          wygrid(3,jcount)=wz
          endif

!1------------- 0 1 0------------------------------
          if(xn.eq.0.0d0.and.yn.ne.0.0d0.and.zn.eq.0.0d0)then
          jcount=jcount+1
          ygrid(1,jcount)=xn
          ygrid(2,jcount)=yn
          ygrid(3,jcount)=zn

          wygrid(1,jcount)=wx
          wygrid(2,jcount)=wy
          wygrid(3,jcount)=wz

!1------------- 0 -1 0------------------------------
          jcount=jcount+1
          ygrid(1,jcount)=xn
          ygrid(2,jcount)=-yn
          ygrid(3,jcount)=zn

          wygrid(1,jcount)=wx
          wygrid(2,jcount)=wy
          wygrid(3,jcount)=wz
          endif

!1------------- 0 0 1------------------------------
          if(xn.eq.0.0d0.and.yn.eq.0.0d0.and.zn.ne.0.0d0)then
          jcount=jcount+1
          ygrid(1,jcount)=xn
          ygrid(2,jcount)=yn
          ygrid(3,jcount)=zn

          wygrid(1,jcount)=wx
          wygrid(2,jcount)=wy
          wygrid(3,jcount)=wz

!1------------- 0 0 -1------------------------------
          jcount=jcount+1
          ygrid(1,jcount)=xn
          ygrid(2,jcount)=yn
          ygrid(3,jcount)=-zn

          wygrid(1,jcount)=wx
          wygrid(2,jcount)=wy
          wygrid(3,jcount)=wz
          endif
!==============================================================

!2------------- 1 1 0------------------------------
          if(xn.ne.0.0d0.and.yn.ne.0.0d0.and.zn.eq.0.0d0)then
          jcount=jcount+1
          ygrid(1,jcount)=xn
          ygrid(2,jcount)=yn
          ygrid(3,jcount)=zn

          wygrid(1,jcount)=wx
          wygrid(2,jcount)=wy
          wygrid(3,jcount)=wz

!2------------- 1 -1 0------------------------------
          jcount=jcount+1
          ygrid(1,jcount)=xn
          ygrid(2,jcount)=-yn
          ygrid(3,jcount)=zn

          wygrid(1,jcount)=wx
          wygrid(2,jcount)=wy
          wygrid(3,jcount)=wz
          
!2------------- -1 1 0------------------------------
          jcount=jcount+1
          ygrid(1,jcount)=-xn
          ygrid(2,jcount)=yn
          ygrid(3,jcount)=zn

          wygrid(1,jcount)=wx
          wygrid(2,jcount)=wy
          wygrid(3,jcount)=wz

!2------------- -1 -1 0------------------------------
          jcount=jcount+1
          ygrid(1,jcount)=-xn
          ygrid(2,jcount)=-yn
          ygrid(3,jcount)=zn

          wygrid(1,jcount)=wx
          wygrid(2,jcount)=wy
          wygrid(3,jcount)=wz

          endif
!==============================================================

!3------------- 1 0 1------------------------------
          if(xn.ne.0.0d0.and.yn.eq.0.0d0.and.zn.ne.0.0d0)then
          jcount=jcount+1
          ygrid(1,jcount)=xn
          ygrid(2,jcount)=yn
          ygrid(3,jcount)=zn

          wygrid(1,jcount)=wx
          wygrid(2,jcount)=wy
          wygrid(3,jcount)=wz

!3------------- 1 0 -1------------------------------
          jcount=jcount+1
          ygrid(1,jcount)=xn
          ygrid(2,jcount)=yn
          ygrid(3,jcount)=-zn

          wygrid(1,jcount)=wx
          wygrid(2,jcount)=wy
          wygrid(3,jcount)=wz
          
!3------------- -1 0 1------------------------------
          jcount=jcount+1
          ygrid(1,jcount)=-xn
          ygrid(2,jcount)=yn
          ygrid(3,jcount)=zn

          wygrid(1,jcount)=wx
          wygrid(2,jcount)=wy
          wygrid(3,jcount)=wz

!3------------- -1 0 -1------------------------------
          jcount=jcount+1
          ygrid(1,jcount)=-xn
          ygrid(2,jcount)=yn
          ygrid(3,jcount)=-zn

          wygrid(1,jcount)=wx
          wygrid(2,jcount)=wy
          wygrid(3,jcount)=wz
          endif
!==============================================================

!4------------- 0 1 1------------------------------
          if(xn.eq.0.0d0.and.yn.ne.0.0d0.and.zn.ne.0.0d0)then
          jcount=jcount+1
          ygrid(1,jcount)=xn
          ygrid(2,jcount)=yn
          ygrid(3,jcount)=zn

          wygrid(1,jcount)=wx
          wygrid(2,jcount)=wy
          wygrid(3,jcount)=wz

!4------------- 0 1 -1------------------------------
          jcount=jcount+1
          ygrid(1,jcount)=xn
          ygrid(2,jcount)=yn
          ygrid(3,jcount)=-zn

          wygrid(1,jcount)=wx
          wygrid(2,jcount)=wy
          wygrid(3,jcount)=wz
          
!4------------- 0 -1 1------------------------------
          jcount=jcount+1
          ygrid(1,jcount)=xn
          ygrid(2,jcount)=-yn
          ygrid(3,jcount)=zn

          wygrid(1,jcount)=wx
          wygrid(2,jcount)=wy
          wygrid(3,jcount)=wz

!4------------- 0 -1 -1------------------------------
          jcount=jcount+1
          ygrid(1,jcount)=xn
          ygrid(2,jcount)=-yn
          ygrid(3,jcount)=-zn

          wygrid(1,jcount)=wx
          wygrid(2,jcount)=wy
          wygrid(3,jcount)=wz
          endif

      enddo

      open(30,file="NGRID",form="formatted")
      rewind(30)
      write(30,*)jcount
      do ii=1,jcount
      write(30,'(6F12.6)')(ygrid(j,ii),j=1,3),(wygrid(k,ii),k=1,3)
      enddo

      print *,'Total GRID points outside WD in first octant:',jcount
      close(30)

      RETURN
      END SUBROUTINE

!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!
!      SUBROUTINE LSGRID
!      INCLUDE 'PARAMS'
!      INCLUDE 'commons.inc'
!      PARAMETER (LM=06)
!      PARAMETER (LMX=3*LM)
!      PARAMETER (MAX_ANG=((LMX+1)*(LMX+2))/2)
!      PARAMETER (MPX=MAX_ANG)
!      PARAMETER (MXPR=MXPOISS)
!
!      DIMENSION SMULTI((LM+1)**2),ILM((LM+1)**2)
!      DIMENSION sump(10),alp(10),rtvec(3)
!      DIMENSION xn1(10),xn2(10),xn3(10),wn1(10),wn2(10),wn3(10)
!
!      COMMON/LATGRID/XVEC(10,MAX_PTS),YVEC(10,MAX_PTS),ZVEC(10,MAX_PTS),
!     &    WXVEC(10,MAX_PTS),WYVEC(10,MAX_PTS),WZVEC(10,MAX_PTS),IGRID
!
!      Character(8) :: Date
!      Character(10) :: Time1,Time2
!      Character(5) :: Zone
!      Integer :: Values(8)
!
!      pi=4.0D0*atan(1.0D0)
!
!!================================================================
!! Gauss quadrature nn-point formula
!! l is the extent of finite sum.
!      nn = 5
!      l  = 10 
!      do m=1,10
!         sump(m)=0.0D0
!         alp(m)=0.0D0
!      enddo
!
!      DO n=1,l-1
!      DO m=2,10
!          sump(m)=sump(m)+1.0D0/(dfloat(n)**dfloat(m))
!      ENDDO
!      ENDDO
!
!      alp(2)=pi*pi/6.0D0-sump(2)
!      alp(2)=1.0D0/alp(2)
!      alp(3)=2.0D0*(1.20205690315959D0-sump(3))
!      alp(3)=1.0D0/alp(3)**(1.0D0/2.0D0)
!      alp(4)=3.0D0*(1.08232323371114D0-sump(4))
!      alp(4)=1.0D0/alp(4)**(1.0D0/3.0D0)
!! Zeta(5) is blowing and for now we approximate zeta(5) with zeta(4)
!      alp(5)=alp(4)
!!      alp(5)=4.0D0*(1.0369277551433699263-sump(5))
!!      alp(5)=1.0D0/alp(5)**(1.0D0/4.0D0)
!      alp(6)=5.0D0*(pi**6.0d0/945.0D0-sump(6))
!      alp(6)=1.0D0/alp(6)**(1.0D0/5.0D0)
!
!      DO ii=7,10
!        alp(ii)=alp(6)
!      ENDDO
!!================================================================
!      igrid=(l**3+3*nn*(l**2)+
!     &      6*(nn**2)*l+6*nn**3)
!
!      if(igrid.GT.MAX_PTS)then
!       PRINT *,'IN LSTRID: INCREASE MAX_PTS',igrid
!       CALL STOPIT
!      endif
!
!      DO jj=1,10
!      DO kk=1,igrid
!        xvec(m,kk)  =0.0d0
!        yvec(m,kk)  =0.0d0
!        zvec(m,kk)  =0.0d0
!        wxvec(m,kk) =0.0d0
!        wyvec(m,kk) =0.0d0
!        wzvec(m,kk) =0.0d0
!      ENDDO
!      ENDDO
!
!      DO ik=1,10
!        xn1(ik)=0.0d0
!        xn2(ik)=0.0d0
!        xn3(ik)=0.0d0
!        wn1(ik)=0.0d0
!        wn2(ik)=0.0d0
!        wn3(ik)=0.0d0
!      ENDDO
!
!      OPEN(299,FILE='SMULTI',FORM='formatted',STATUS='OLD')
!      REWIND(299)
!      READ(299,*)ILMTOT
!      m_old=0
!
!      DO 505 kk=1,ILMTOT
!      READ(299,*)ILM(kk),SMULTI(kk)
!
!      mm=ILM(kk) 
!      m=ceiling(dsqrt(dfloat(mm)))
!      if(m.LE.m_old)goto 505
!
!      iicount=0
!
!      iwx=1
!      iwy=1
!      iwz=1
!      if(NXYZ.ge.1)iwz=l
!      if(NXYZ.ge.2)iwy=l
!      if(NXYZ.ge.3)iwx=l
!
!      do n1=0,iwx-1
!      do n2=0,iwy-1
!      do n3=0,iwz-1
!         if(n1.eq.0.and.n2.eq.0.and.n3.eq.0)goto 800
!         iicount=iicount+1
!         wxvec(m,iicount)=1.0d0
!         wyvec(m,iicount)=1.0d0
!         wzvec(m,iicount)=1.0d0
!         xvec(m,iicount)=dfloat(n1)
!         yvec(m,iicount)=dfloat(n2)
!         zvec(m,iicount)=dfloat(n3)
!  800    continue
!
!         enddo ! n3 loop ends
!!         if(nxyz.le.2)goto 504
!!==========================================================================
!
!            dlimn30 = 0.0D0
!            dlimn31 = 1.0D0/alp(m)
!
!            call gaussp(dlimn30,dlimn31,nn,xn3,wn3)
!            do k=1,nn
!
!               iicount=iicount+1
!
!               wxvec(m,iicount)=1.0d0
!               wyvec(m,iicount)=1.0d0
!               wzvec(m,iicount)=wn3(k)/(xn3(k)**2.0d0)
!
!               xvec(m,iicount)=dfloat(n1)
!               yvec(m,iicount)=dfloat(n2)
!               zvec(m,iicount)=1.0d0/xn3(k)
!            enddo
!
!  504   continue
!        enddo !n2 loop ends      
!!------------------------------------------------------------------------            
!        if(nxyz.le.1)goto 503
!
!        dlimn20=0.0D0
!        dlimn21=pi/4.0D0
!        call gaussp(dlimn20,dlimn21,nn,xn2,wn2)
!        do j=1,nn
!           dlimn30=0.0D0
!           dlimn31=1.0D0/(alp(m)*cos(xn2(j)))
!           call gaussp(dlimn30,dlimn31,nn,xn3,wn3)
!           do k=1,nn
!
!              rr=xn3(k)
!              vv=rr*cos(xn2(j))
!              ww=rr*sin(xn2(j))
!
!              iicount=iicount+1
!
!              wxvec(m,iicount)=1.0d0
!              wyvec(m,iicount)=wn2(j)*(dsqrt(rr)/vv**2.0d0)
!              wzvec(m,iicount)=wn3(k)*(dsqrt(rr)/ww**2.0d0)
!              xvec(m,iicount)=dfloat(n1)
!              yvec(m,iicount)=1.0d0/vv
!              zvec(m,iicount)=1.0d0/ww
!
!!========================================================================              
!!             Term 2
!              vv=rr*sin(xn2(j))
!              ww=rr*cos(xn2(j))
!
!              iicount=iicount+1
!
!              wxvec(m,iicount)=1.0d0
!              wyvec(m,iicount)=wn2(j)*(dsqrt(rr)/vv**2.0d0)
!              wzvec(m,iicount)=wn3(k)*(dsqrt(rr)/ww**2.0d0)
!              xvec(m,iicount)=dfloat(n1)
!              yvec(m,iicount)=1.0d0/vv
!              zvec(m,iicount)=1.0d0/ww
!           enddo
!        enddo
! 503  continue
!      enddo ! n1 loop ends.       
!!=====================================================================              
!      do n3=0,l-1
!        do n1=0,l-1
!          dlimn20 = 0.0D0
!          dlimn21 = 1.0D0/alp(m)
!          call gaussp(dlimn20,dlimn21,nn,xn2,wn2)
!          do j=1,nn
!            iicount=iicount+1
!
!            wxvec(m,iicount)=1.0d0
!            wyvec(m,iicount)=wn2(j)/(xn2(j)**2.0d0)
!            wzvec(m,iicount)=1.0d0
!            xvec(m,iicount)=dfloat(n1)
!            yvec(m,iicount)=1.0d0/xn2(j)
!            zvec(m,iicount)=dfloat(n3)
!          enddo
!         enddo ! n1 loop ends
!!=======================================================================
!
!         dlimn20 = 0.0D0
!         dlimn21 = pi/4.0D0
!         call gaussp(dlimn20,dlimn21,nn,xn2,wn2)
!         do j=1,nn
!            dlimn10 = 0.0D0
!            dlimn11 = 1.0D0/alp(m)*cos(xn2(j))
!            call gaussp(dlimn10,dlimn11,nn,xn1,wn1)
!            do i=1,nn
!              rr=xn1(i)
!              uu=rr*cos(xn2(j))
!              vv=rr*sin(xn2(j))
!              iicount=iicount+1
!
!              wxvec(m,iicount)=wn1(i)*(dsqrt(rr)/uu**2.0d0)
!              wyvec(m,iicount)=wn2(j)*(dsqrt(rr)/vv**2.0d0)
!              wzvec(m,iicount)=1.0d0
!              xvec(m,iicount)=1.0d0/uu
!              yvec(m,iicount)=1.0d0/vv
!              zvec(m,iicount)=dfloat(n3)
!
!!=====================================================================
!              ! PART 2
!
!              uu=rr*sin(xn2(j))
!              vv=rr*cos(xn2(j))
!
!              iicount=iicount+1
!
!              wxvec(m,iicount)=wn1(i)*(dsqrt(rr)/uu**2.0d0)
!              wyvec(m,iicount)=wn2(j)*(dsqrt(rr)/vv**2.0d0)
!              wzvec(m,iicount)=1.0d0
!              xvec(m,iicount)=1.0d0/uu
!              yvec(m,iicount)=1.0d0/vv
!              zvec(m,iicount)=dfloat(n3)
!
!            enddo
!         enddo
!        enddo ! n3 loop      
!!========================================================================
!
!        do n2=0,l-1
!           do n3=0,l-1
!              dlimn10 = 0.0D0
!              dlimn11 = 1.0D0/alp(m)
!              call gaussp(dlimn10,dlimn11,nn,xn1,wn1)
!              do i=1,nn
!                 iicount=iicount+1
!
!                 wxvec(m,iicount)=wn1(i)/(xn1(i)**2.0d0)
!                 wyvec(m,iicount)=1.0d0
!                 wzvec(m,iicount)=1.0d0
!                 xvec(m,iicount)=1.0d0/xn1(i)
!                 yvec(m,iicount)=dfloat(n2)
!                 zvec(m,iicount)=dfloat(n3)
!               enddo
!           enddo ! n3 loop ends
!
!!====================================================================================
!
!           dlimn10 = 0.0D0
!           dlimn11 = pi/4.0D0
!           call gaussp(dlimn10,dlimn11,nn,xn1,wn1)
!           do i=1,nn
!              dlimn30 = 0.0D0
!              dlimn31 = 1.0D0/(alp(m)*cos(xn1(i)))
!              call gaussp(dlimn30,dlimn31,nn,xn3,wn3)
!              do k=1,nn
!                 rr=xn3(k)
!                 uu=rr*cos(xn1(i))
!                 ww=rr*sin(xn1(i))
!
!                 iicount=iicount+1
!
!                 wxvec(m,iicount)=wn1(i)*(dsqrt(rr)/uu**2.0d0)
!                 wyvec(m,iicount)=1.0d0
!                 wzvec(m,iicount)=wn3(k)*(dsqrt(rr)/ww**2.0d0)
!                 xvec(m,iicount)=1.0d0/uu
!                 yvec(m,iicount)=dfloat(n2)
!                 zvec(m,iicount)=1.0d0/ww
!
!!=========================================================================
!                 uu=rr*sin(xn1(i))
!                 ww=rr*cos(xn1(i))
!
!                 iicount=iicount+1
!
!                 wxvec(m,iicount)=wn1(i)*(dsqrt(rr)/uu**2.0d0)
!                 wyvec(m,iicount)=1.0d0
!                 wzvec(m,iicount)=wn3(k)*(dsqrt(rr)/ww**2.0d0)
!                 xvec(m,iicount)=1.0d0/uu
!                 yvec(m,iicount)=dfloat(n2)
!                 zvec(m,iicount)=1.0d0/ww
!               enddo
!            enddo
!        enddo !n2 loop
!
!!======================================================================
!        ! Triple Integration
!!======================================================================
!        dlimn30 = 0.0D0
!        dlimn31 = pi/4.0D0
!        call gaussp(dlimn30,dlimn31,nn,xn3,wn3)
!        do k=1,nn
!           dlimn20 = atan(1.0D0/sin(xn3(k)))
!           dlimn21 = pi/2.0D0
!           call gaussp(dlimn20,dlimn21,nn,xn2,wn2)
!           do j=1,nn
!              dlimn10 = 0.0D0
!              dlimn11 = 1.0D0/(alp(m)*sin(xn2(j))*cos(xn3(k)))
!              call gaussp(dlimn10,dlimn11,nn,xn1,wn1)
!              do i=1,nn
!! Term 1
!                 rr=xn1(i)
!                 wttfact=(rr*rr*sin(xn2(j)))**(1.0d0/3.0d0)
!  
!                 uu=rr*sin(xn2(j))*cos(xn3(k))
!                 vv=rr*sin(xn2(j))*sin(xn3(k))
!                 ww=rr*cos(xn2(j))
!  
!                 iicount=iicount+1
!
!                 wxvec(m,iicount)=(wn1(i)*wttfact)/(uu**2.0d0)
!                 wyvec(m,iicount)=(wn2(j)*wttfact)/(vv**2.0d0)
!                 wzvec(m,iicount)=(wn3(k)*wttfact)/(ww**2.0d0)
!                 xvec(m,iicount)=1.0d0/uu
!                 yvec(m,iicount)=1.0d0/vv
!                 zvec(m,iicount)=1.0d0/ww
!! Term 2
!                 uu=rr*sin(xn2(j))*sin(xn3(k))
!                 vv=rr*sin(xn2(j))*cos(xn3(k))
!                 ww=rr*cos(xn2(j))
!  
!                 iicount=iicount+1
!
!                 wxvec(m,iicount)=(wn1(i)*wttfact)/(uu**2.0d0)
!                 wyvec(m,iicount)=(wn2(j)*wttfact)/(vv**2.0d0)
!                 wzvec(m,iicount)=(wn3(k)*wttfact)/(ww**2.0d0)
!                 xvec(m,iicount)=1.0d0/uu
!                 yvec(m,iicount)=1.0d0/vv
!                 zvec(m,iicount)=1.0d0/ww
!! Term   3
!                 uu=rr*sin(xn2(j))*sin(xn3(k))
!                 vv=rr*cos(xn2(j))
!                 ww=rr*sin(xn2(j))*cos(xn3(k))
!   
!                 iicount=iicount+1
!
!                 wxvec(m,iicount)=(wn1(i)*wttfact)/(uu**2.0d0)
!                 wyvec(m,iicount)=(wn2(j)*wttfact)/(vv**2.0d0)
!                 wzvec(m,iicount)=(wn3(k)*wttfact)/(ww**2.0d0)
!                 xvec(m,iicount)=1.0d0/uu
!                 yvec(m,iicount)=1.0d0/vv
!                 zvec(m,iicount)=1.0d0/ww
!! Term 4
!                 uu=rr*sin(xn2(j))*cos(xn3(k))
!                 vv=rr*cos(xn2(j))
!                 ww=rr*sin(xn2(j))*sin(xn3(k))
!
!                 iicount=iicount+1
!
!                 wxvec(m,iicount)=(wn1(i)*wttfact)/(uu**2.0d0)
!                 wyvec(m,iicount)=(wn2(j)*wttfact)/(vv**2.0d0)
!                 wzvec(m,iicount)=(wn3(k)*wttfact)/(ww**2.0d0)
!                 xvec(m,iicount)=1.0d0/uu
!                 yvec(m,iicount)=1.0d0/vv
!                 zvec(m,iicount)=1.0d0/ww
!! Term 5
!                 xx=rr*cos(xn2(j))
!                 vv=rr*sin(xn2(j))*cos(xn3(k))
!                 ww=rr*sin(xn2(j))*sin(xn3(k))
!
!                 iicount=iicount+1
!
!                 wxvec(m,iicount)=(wn1(i)*wttfact)/(uu**2.0d0)
!                 wyvec(m,iicount)=(wn2(j)*wttfact)/(vv**2.0d0)
!                 wzvec(m,iicount)=(wn3(k)*wttfact)/(ww**2.0d0)
!                 xvec(m,iicount)=1.0d0/uu
!                 yvec(m,iicount)=1.0d0/vv
!                 zvec(m,iicount)=1.0d0/ww
!! Term 6
!                 xx=rr*cos(xn2(j))
!                 vv=rr*sin(xn2(j))*sin(xn3(k))
!                 ww=rr*sin(xn2(j))*cos(xn3(k))
!   
!                 iicount=iicount+1
!
!                 wxvec(m,iicount)=(wn1(i)*wttfact)/(uu**2.0d0)
!                 wyvec(m,iicount)=(wn2(j)*wttfact)/(vv**2.0d0)
!                 wzvec(m,iicount)=(wn3(k)*wttfact)/(ww**2.0d0)
!                 xvec(m,iicount)=1.0d0/uu
!                 yvec(m,iicount)=1.0d0/vv
!                 zvec(m,iicount)=1.0d0/ww
!          enddo
!        enddo       
!       enddo       
!
!      m_old=m
!
!  505 CONTINUE ! m loop started at the top ends here 
!      CLOSE(299)
!
!!      do ik=1,5
!!      do jc=1,iicount
!!       print *,'x y z', xvec(ik,jc),yvec(ik,jc),zvec(ik,jc)
!!      enddo
!!      enddo
!
!       if(iicount.gt.igrid)then
!        print *,'iicount.GT.igrid',iicount,igrid
!        CALL STOPIT
!       endif
!       print *,''
!       igrid=iicount
!       print *,'Total GRID points outside WD in first octant:',igrid
!
!       RETURN
!      END SUBROUTINE LSGRID
!
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      SUBROUTINE LAT_SUM(r0,sumdirect)
      INCLUDE 'PARAMS'
      INCLUDE 'commons.inc'
      PARAMETER (LM=06)
      PARAMETER (LMX=3*LM)
      PARAMETER (MAX_ANG=((LMX+1)*(LMX+2))/2)
      PARAMETER (MPX=MAX_ANG)
      PARAMETER (MXPR=MXPOISS)
      PARAMETER (l=10)
      PARAMETER (nn=5)
      PARAMETER (igrid=(l**2+2*nn*l+2*nn**2))
      COMMON/LATGRID/XVEC(10,MAX_PTS),YVEC(10,MAX_PTS),ZVEC(10,MAX_PTS),
     &    WXVEC(10,MAX_PTS),WYVEC(10,MAX_PTS),WZVEC(10,MAX_PTS)
      DIMENSION SMULTI((LM+1)**2),ILM((LM+1)**2)
      DIMENSION rtvec(3),TVEC(3,3),r0(3),ccenter(3)
      DIMENSION unitvec2(3,8*igrid),POLY2(8*igrid,(LM+1)**2,LM+1)
      DIMENSION unitvec(3,8*igrid,LM+1),denom(8*igrid,LM+1),
     &      wtt(8*igrid,LM+1),POLY(MPX,(LM+1)**2)

      Character(8) :: Date
      Character(10) :: Time1,Time2
      Character(5) :: Zone
      Integer :: Values(8)

      pi = 4.0D0*atan(1.0D0)

      ccenter=0.0D0
      r0(1)=r0(1)-ccenter(1)
      r0(2)=r0(2)-ccenter(2)
      r0(3)=r0(3)-ccenter(3)

      unitvec=0.0D0
      m_old=0

!      OPEN(299,FILE='SMULTI',FORM='formatted',STATUS='OLD')
!      REWIND(299)
!      !READ(299,*)(ccenter(ii),ii=1,3)
!      READ(299,*)ILMTOT
!
!      DO 505 kk=1,ILMTOT
!      READ(299,*)ILM(kk),SMULTI(kk)
!
!      mm=ILM(kk) 
!      m=ceiling(dsqrt(dfloat(mm)))
!
!      if(m.LE.m_old) goto 505

      open(30,file="NGRID")
      read(30,*)jcount

      DO ii=1,jcount
          read(30,*)xnn,ynn,znn,wx,wy,wz

          wtt(ii,m)=wx*wy*wz
!          denom(ii,m)=dsqrt((r0(1)-xnn)**2+(r0(2)-ynn)**2+(r0(3)-znn)**2)
          do ix=1,3
             unitvec(ix,jjcount,m)=rtvec(ix)/dsqrt(denom(jjcount,m))
          enddo

          !print 507,jjcount,xnn2,ynn2,znn2,absrtvec,wtt(jjcount,m)

!1------------- -1 0 0------------------------------
          xnn2=-xnn
          ynn2=ynn
          znn2=znn

          jjcount=jjcount+1
          wtt(jjcount,m)=wxvec(m,ii)*wyvec(m,ii)*wzvec(m,ii)
          rtvec(:)=(r0(:)-xnn2*TVEC(:,1)-ynn2*TVEC(:,2)-znn2*TVEC(:,3))
          denom(jjcount,m)=(sum(rtvec(:)*rtvec(:)))!**(dfloat(m)/2.0d0)
          do ix=1,3
             unitvec(ix,jjcount,m)=rtvec(ix)/dsqrt(denom(jjcount,m))
          enddo
          !print 507,jjcount,xnn2,ynn2,znn2,absrtvec,wtt(jjcount,m)

!1------------- 0 1 0------------------------------
        if(xnn.eq.0.0d0.and.ynn.ne.0.0d0.and.znn.eq.0.0d0)then

          xnn2=xnn
          ynn2=ynn
          znn2=znn

          jjcount=jjcount+1
          wtt(jjcount,m)=wxvec(m,ii)*wyvec(m,ii)*wzvec(m,ii)
          rtvec(:)=(r0(:)-xnn2*TVEC(:,1)-ynn2*TVEC(:,2)-znn2*TVEC(:,3))
          denom(jjcount,m)=(sum(rtvec(:)*rtvec(:)))!**(dfloat(m)/2.0d0)
          do ix=1,3
             unitvec(ix,jjcount,m)=rtvec(ix)/dsqrt(denom(jjcount,m))
          enddo
          !print 507,jjcount,xnn2,ynn2,znn2,absrtvec,wtt(jjcount,m)

!1------------- 0 -1 0------------------------------
          xnn2=xnn
          ynn2=-ynn
          znn2=znn

          jjcount=jjcount+1
          wtt(jjcount,m)=wxvec(m,ii)*wyvec(m,ii)*wzvec(m,ii)
          rtvec(:)=(r0(:)-xnn2*TVEC(:,1)-ynn2*TVEC(:,2)-znn2*TVEC(:,3))
          denom(jjcount,m)=(sum(rtvec(:)*rtvec(:)))!**(dfloat(m)/2.0d0)
          do ix=1,3
             unitvec(ix,jjcount,m)=rtvec(ix)/dsqrt(denom(jjcount,m))
          enddo
          !print 507,jjcount,xnn2,ynn2,znn2,absrtvec,wtt(jjcount,m)
        endif

!1------------- 0 0 1------------------------------
        if(xnn.eq.0.0d0.and.ynn.eq.0.0d0.and.znn.ne.0.0d0)then

          xnn2=xnn
          ynn2=ynn
          znn2=znn

          jjcount=jjcount+1
          wtt(jjcount,m)=wxvec(m,ii)*wyvec(m,ii)*wzvec(m,ii)
          rtvec(:)=(r0(:)-xnn2*TVEC(:,1)-ynn2*TVEC(:,2)-znn2*TVEC(:,3))
          denom(jjcount,m)=(sum(rtvec(:)*rtvec(:)))!**(dfloat(m)/2.0d0)
          do ix=1,3
             unitvec(ix,jjcount,m)=rtvec(ix)/dsqrt(denom(jjcount,m))
          enddo
          !print 507,jjcount,xnn2,ynn2,znn2,absrtvec,wtt(jjcount,m)

!1------------- 0 0 -1------------------------------
          xnn2=xnn
          ynn2=ynn
          znn2=-znn

          jjcount=jjcount+1
          wtt(jjcount,m)=wxvec(m,ii)*wyvec(m,ii)*wzvec(m,ii)
          rtvec(:)=(r0(:)-xnn2*TVEC(:,1)-ynn2*TVEC(:,2)-znn2*TVEC(:,3))
          denom(jjcount,m)=(sum(rtvec(:)*rtvec(:)))!**(dfloat(m)/2.0d0)
          do ix=1,3
             unitvec(ix,jjcount,m)=rtvec(ix)/dsqrt(denom(jjcount,m))
          enddo
          !print 507,jjcount,xnn2,ynn2,znn2,absrtvec,wtt(jjcount,m)
        endif

!1------------- 0 1 1------------------------------
        if(xnn.eq.0.0d0.and.ynn.ne.0.0d0.and.znn.ne.0.0d0)then

          xnn2=xnn
          ynn2=ynn
          znn2=znn

          jjcount=jjcount+1
          wtt(jjcount,m)=wxvec(m,ii)*wyvec(m,ii)*wzvec(m,ii)
          rtvec(:)=(r0(:)-xnn2*TVEC(:,1)-ynn2*TVEC(:,2)-znn2*TVEC(:,3))
          denom(jjcount,m)=(sum(rtvec(:)*rtvec(:)))!**(dfloat(m)/2.0d0)
          do ix=1,3
             unitvec(ix,jjcount,m)=rtvec(ix)/dsqrt(denom(jjcount,m))
          enddo
          !print 507,jjcount,xnn2,ynn2,znn2,absrtvec,wtt(jjcount,m)

!1------------- 0 -1 1------------------------------
          xnn2=xnn
          ynn2=-ynn
          znn2=znn

          jjcount=jjcount+1
          wtt(jjcount,m)=wxvec(m,ii)*wyvec(m,ii)*wzvec(m,ii)
          rtvec(:)=(r0(:)-xnn2*TVEC(:,1)-ynn2*TVEC(:,2)-znn2*TVEC(:,3))
          denom(jjcount,m)=(sum(rtvec(:)*rtvec(:)))!**(dfloat(m)/2.0d0)
          do ix=1,3
             unitvec(ix,jjcount,m)=rtvec(ix)/dsqrt(denom(jjcount,m))
          enddo
          !print 507,jjcount,xnn2,ynn2,znn2,absrtvec,wtt(jjcount,m)

!1------------- 0 1 -1------------------------------
          xnn2=xnn
          ynn2=ynn
          znn2=-znn

          jjcount=jjcount+1
          wtt(jjcount,m)=wxvec(m,ii)*wyvec(m,ii)*wzvec(m,ii)
          rtvec(:)=(r0(:)-xnn2*TVEC(:,1)-ynn2*TVEC(:,2)-znn2*TVEC(:,3))
          denom(jjcount,m)=(sum(rtvec(:)*rtvec(:)))!**(dfloat(m)/2.0d0)
          do ix=1,3
             unitvec(ix,jjcount,m)=rtvec(ix)/dsqrt(denom(jjcount,m))
          enddo
          !print 507,jjcount,xnn2,ynn2,znn2,absrtvec,wtt(jjcount,m)

!1------------- 0 -1 -1------------------------------
          xnn2=xnn
          ynn2=-ynn
          znn2=-znn

          jjcount=jjcount+1
          wtt(jjcount,m)=wxvec(m,ii)*wyvec(m,ii)*wzvec(m,ii)
          rtvec(:)=(r0(:)-xnn2*TVEC(:,1)-ynn2*TVEC(:,2)-znn2*TVEC(:,3))
          denom(jjcount,m)=(sum(rtvec(:)*rtvec(:)))!**(dfloat(m)/2.0d0)
          do ix=1,3
             unitvec(ix,jjcount,m)=rtvec(ix)/dsqrt(denom(jjcount,m))
          enddo
          !print 507,jjcount,xnn2,ynn2,znn2,absrtvec,wtt(jjcount,m)

        endif
!1------------- 1 0 1------------------------------
        if(xnn.ne.0.0d0.and.ynn.eq.0.0d0.and.znn.ne.0.0d0)then

          xnn2=xnn
          ynn2=ynn
          znn2=znn

          jjcount=jjcount+1
          wtt(jjcount,m)=wxvec(m,ii)*wyvec(m,ii)*wzvec(m,ii)
          rtvec(:)=(r0(:)-xnn2*TVEC(:,1)-ynn2*TVEC(:,2)-znn2*TVEC(:,3))
          denom(jjcount,m)=(sum(rtvec(:)*rtvec(:)))!**(dfloat(m)/2.0d0)
          do ix=1,3
             unitvec(ix,jjcount,m)=rtvec(ix)/dsqrt(denom(jjcount,m))
          enddo
          !print 507,jjcount,xnn2,ynn2,znn2,absrtvec,wtt(jjcount,m)

!1------------- -1 0 1------------------------------
          xnn2=-xnn
          ynn2=ynn
          znn2=znn

          jjcount=jjcount+1
          wtt(jjcount,m)=wxvec(m,ii)*wyvec(m,ii)*wzvec(m,ii)
          rtvec(:)=(r0(:)-xnn2*TVEC(:,1)-ynn2*TVEC(:,2)-znn2*TVEC(:,3))
          denom(jjcount,m)=(sum(rtvec(:)*rtvec(:)))!**(dfloat(m)/2.0d0)
          do ix=1,3
             unitvec(ix,jjcount,m)=rtvec(ix)/dsqrt(denom(jjcount,m))
          enddo
          !print 507,jjcount,xnn2,ynn2,znn2,absrtvec,wtt(jjcount,m)

!1------------- 1 0 -1------------------------------
          xnn2=xnn
          ynn2=ynn
          znn2=-znn

          jjcount=jjcount+1
          wtt(jjcount,m)=wxvec(m,ii)*wyvec(m,ii)*wzvec(m,ii)
          rtvec(:)=(r0(:)-xnn2*TVEC(:,1)-ynn2*TVEC(:,2)-znn2*TVEC(:,3))
          denom(jjcount,m)=(sum(rtvec(:)*rtvec(:)))!**(dfloat(m)/2.0d0)
          do ix=1,3
             unitvec(ix,jjcount,m)=rtvec(ix)/dsqrt(denom(jjcount,m))
          enddo
          !print 507,jjcount,xnn2,ynn2,znn2,absrtvec,wtt(jjcount,m)

!1------------- -1 0 -1------------------------------
          xnn2=-xnn
          ynn2=ynn
          znn2=-znn

          jjcount=jjcount+1
          wtt(jjcount,m)=wxvec(m,ii)*wyvec(m,ii)*wzvec(m,ii)
          rtvec(:)=(r0(:)-xnn2*TVEC(:,1)-ynn2*TVEC(:,2)-znn2*TVEC(:,3))
          denom(jjcount,m)=(sum(rtvec(:)*rtvec(:)))!**(dfloat(m)/2.0d0)
          do ix=1,3
             unitvec(ix,jjcount,m)=rtvec(ix)/dsqrt(denom(jjcount,m))
          enddo
          !print 507,jjcount,xnn2,ynn2,znn2,absrtvec,wtt(jjcount,m)

        endif
!1------------- 1 1 0------------------------------
        if(xnn.ne.0.0d0.and.ynn.ne.0.0d0.and.znn.eq.0.0d0)then

          xnn2=xnn
          ynn2=ynn
          znn2=znn

          jjcount=jjcount+1
          wtt(jjcount,m)=wxvec(m,ii)*wyvec(m,ii)*wzvec(m,ii)
          rtvec(:)=(r0(:)-xnn2*TVEC(:,1)-ynn2*TVEC(:,2)-znn2*TVEC(:,3))
          denom(jjcount,m)=(sum(rtvec(:)*rtvec(:)))!**(dfloat(m)/2.0d0)
          do ix=1,3
             unitvec(ix,jjcount,m)=rtvec(ix)/dsqrt(denom(jjcount,m))
          enddo
          !print 507,jjcount,xnn2,ynn2,znn2,absrtvec,wtt(jjcount,m)

!1------------- -1 1 0------------------------------
          xnn2=-xnn
          ynn2=ynn
          znn2=znn

          jjcount=jjcount+1
          wtt(jjcount,m)=wxvec(m,ii)*wyvec(m,ii)*wzvec(m,ii)
          rtvec(:)=(r0(:)-xnn2*TVEC(:,1)-ynn2*TVEC(:,2)-znn2*TVEC(:,3))
          denom(jjcount,m)=(sum(rtvec(:)*rtvec(:)))!**(dfloat(m)/2.0d0)
          do ix=1,3
             unitvec(ix,jjcount,m)=rtvec(ix)/dsqrt(denom(jjcount,m))
          enddo
          !print 507,jjcount,xnn2,ynn2,znn2,absrtvec,wtt(jjcount,m)

!1------------- 1 -1 0------------------------------
          xnn2=xnn
          ynn2=-ynn
          znn2=znn

          jjcount=jjcount+1
          wtt(jjcount,m)=wxvec(m,ii)*wyvec(m,ii)*wzvec(m,ii)
          rtvec(:)=(r0(:)-xnn2*TVEC(:,1)-ynn2*TVEC(:,2)-znn2*TVEC(:,3))
          denom(jjcount,m)=(sum(rtvec(:)*rtvec(:)))!**(dfloat(m)/2.0d0)
          do ix=1,3
             unitvec(ix,jjcount,m)=rtvec(ix)/dsqrt(denom(jjcount,m))
          enddo
          !print 507,jjcount,xnn2,ynn2,znn2,absrtvec,wtt(jjcount,m)

!1------------- -1 -1 0------------------------------
          xnn2=-xnn
          ynn2=-ynn
          znn2=znn

          jjcount=jjcount+1
          wtt(jjcount,m)=wxvec(m,ii)*wyvec(m,ii)*wzvec(m,ii)
          rtvec(:)=(r0(:)-xnn2*TVEC(:,1)-ynn2*TVEC(:,2)-znn2*TVEC(:,3))
          denom(jjcount,m)=(sum(rtvec(:)*rtvec(:)))!**(dfloat(m)/2.0d0)
          do ix=1,3
             unitvec(ix,jjcount,m)=rtvec(ix)/dsqrt(denom(jjcount,m))
          enddo
          !print 507,jjcount,xnn2,ynn2,znn2,absrtvec,wtt(jjcount,m)
        endif

        if(xnn.ne.0.0d0.and.ynn.ne.0.0d0.and.znn.ne.0.0d0)then
!1------------- 1 1 1------------------------------
          xnn2=xnn
          ynn2=ynn
          znn2=znn

          jjcount=jjcount+1
          wtt(jjcount,m)=wxvec(m,ii)*wyvec(m,ii)*wzvec(m,ii)
          rtvec(:)=(r0(:)-xnn2*TVEC(:,1)-ynn2*TVEC(:,2)-znn2*TVEC(:,3))
          denom(jjcount,m)=(sum(rtvec(:)*rtvec(:)))!**(dfloat(m)/2.0d0)
          do ix=1,3
             unitvec(ix,jjcount,m)=rtvec(ix)/dsqrt(denom(jjcount,m))
          enddo

          !print 507,jjcount,xnn2,ynn2,znn2,absrtvec,wtt(jjcount,m)
!2------------- 1 1 -1------------------------------
          xnn2=xnn
          ynn2=ynn
          znn2=-znn

          jjcount=jjcount+1
          wtt(jjcount,m)=wxvec(m,ii)*wyvec(m,ii)*wzvec(m,ii)
          rtvec(:)=(r0(:)-xnn2*TVEC(:,1)-ynn2*TVEC(:,2)-znn2*TVEC(:,3))
          denom(jjcount,m)=(sum(rtvec(:)*rtvec(:)))!**(dfloat(m)/2.0d0)
          do ix=1,3
             unitvec(ix,jjcount,m)=rtvec(ix)/dsqrt(denom(jjcount,m))
          enddo
          !print 507,jjcount,xnn2,ynn2,znn2,absrtvec,wtt(jjcount,m)

!3------------- 1 -1 1------------------------------
          xnn2=xnn
          ynn2=-ynn
          znn2=znn

          jjcount=jjcount+1
          wtt(jjcount,m)=wxvec(m,ii)*wyvec(m,ii)*wzvec(m,ii)
          rtvec(:)=(r0(:)-xnn2*TVEC(:,1)-ynn2*TVEC(:,2)-znn2*TVEC(:,3))
          denom(jjcount,m)=(sum(rtvec(:)*rtvec(:)))!**(dfloat(m)/2.0d0)
          do ix=1,3
             unitvec(ix,jjcount,m)=rtvec(ix)/dsqrt(denom(jjcount,m))
          enddo
          !print 507,jjcount,xnn2,ynn2,znn2,absrtvec,wtt(jjcount,m)

!4------------- -1 1 1------------------------------
          xnn2=-xnn
          ynn2=ynn
          znn2=znn

          jjcount=jjcount+1
          wtt(jjcount,m)=wxvec(m,ii)*wyvec(m,ii)*wzvec(m,ii)
          rtvec(:)=(r0(:)-xnn2*TVEC(:,1)-ynn2*TVEC(:,2)-znn2*TVEC(:,3))
          denom(jjcount,m)=(sum(rtvec(:)*rtvec(:)))!**(dfloat(m)/2.0d0)
          do ix=1,3
             unitvec(ix,jjcount,m)=rtvec(ix)/dsqrt(denom(jjcount,m))
          enddo
          !print 507,jjcount,xnn2,ynn2,znn2,absrtvec,wtt(jjcount,m)

!5------------- -1 -1 1------------------------------
          xnn2=-xnn
          ynn2=-ynn
          znn2=znn

          jjcount=jjcount+1
          wtt(jjcount,m)=wxvec(m,ii)*wyvec(m,ii)*wzvec(m,ii)
          rtvec(:)=(r0(:)-xnn2*TVEC(:,1)-ynn2*TVEC(:,2)-znn2*TVEC(:,3))
          denom(jjcount,m)=(sum(rtvec(:)*rtvec(:)))!**(dfloat(m)/2.0d0)
          do ix=1,3
             unitvec(ix,jjcount,m)=rtvec(ix)/dsqrt(denom(jjcount,m))
          enddo
          !print 507,jjcount,xnn2,ynn2,znn2,absrtvec,wtt(jjcount,m)

!6------------- -1 1 -1------------------------------
          xnn2=-xnn
          ynn2=ynn
          znn2=-znn

          jjcount=jjcount+1
          wtt(jjcount,m)=wxvec(m,ii)*wyvec(m,ii)*wzvec(m,ii)
          rtvec(:)=(r0(:)-xnn2*TVEC(:,1)-ynn2*TVEC(:,2)-znn2*TVEC(:,3))
          denom(jjcount,m)=(sum(rtvec(:)*rtvec(:)))!**(dfloat(m)/2.0d0)
          do ix=1,3
             unitvec(ix,jjcount,m)=rtvec(ix)/dsqrt(denom(jjcount,m))
          enddo
          !print 507,jjcount,xnn2,ynn2,znn2,absrtvec,wtt(jjcount,m)

!7-------------1 -1 -1------------------------------

          xnn2=xnn
          ynn2=-ynn
          znn2=-znn

          jjcount=jjcount+1
          wtt(jjcount,m)=wxvec(m,ii)*wyvec(m,ii)*wzvec(m,ii)
          rtvec(:)=(r0(:)-xnn2*TVEC(:,1)-ynn2*TVEC(:,2)-znn2*TVEC(:,3))
          denom(jjcount,m)=(sum(rtvec(:)*rtvec(:)))!**(dfloat(m)/2.0d0)
          do ix=1,3
             unitvec(ix,jjcount,m)=rtvec(ix)/dsqrt(denom(jjcount,m))
          enddo
          !print 507,jjcount,xnn2,ynn2,znn2,absrtvec,wtt(jjcount,m)

!8------------- -1 -1 -1------------------------------
          xnn2=-xnn
          ynn2=-ynn
          znn2=-znn

          jjcount=jjcount+1
          wtt(jjcount,m)=wxvec(m,ii)*wyvec(m,ii)*wzvec(m,ii)
          rtvec(:)=(r0(:)-xnn2*TVEC(:,1)-ynn2*TVEC(:,2)-znn2*TVEC(:,3))
          denom(jjcount,m)=(sum(rtvec(:)*rtvec(:)))!**(dfloat(m)/2.0d0)
          do ix=1,3
             unitvec(ix,jjcount,m)=rtvec(ix)/dsqrt(denom(jjcount,m))
          enddo
          !print 507,jjcount,xnn2,ynn2,znn2,absrtvec,wtt(jjcount,m)
        endif
  151   CONTINUE
        ENDDO
        m_old=m
  505 CONTINUE

      CLOSE(299)
! REARRAGE jjcount in increasing order of radius

!      do jj=1,jjcount
!        !print *,(unitvec(ix,jj,m),ix=1,3)
!      enddo

      !PRINT *,'GRID SIZE for each m (l+1) value:',jjcount

      POLY=0.0d0
      POLY2=0.0d0
      m_old=0
      DO 501 kk=1,ILMTOT
       mm=ILM(kk)
       m=ceiling(dsqrt(dfloat(mm)))
       if(m.LE.m_old)goto 501
       do ix=1,3
       unitvec2(ix,jjcount)=unitvec(ix,jjcount,m)
       enddo
       CALL HARMONICS(MPX,jjcount,LM,unitvec2,POLY,NPOLY)

       do jj=1,jjcount
       do ii=1,NPOLY
        POLY2(jj,ii,m)=POLY(jj,ii)
       enddo
       enddo
       m_old=m
  501 CONTINUE

      sumdirect=0.0d0

      do jj=1,jjcount
      DO 504 kk=1,ILMTOT
        mm=ILM(kk) 
        m=ceiling(dsqrt(dfloat(lm)))
        sumdirect=sumdirect+(SMULTI(kk)*POLY2(jj,mm,m)*wtt(jj,m))
     &             /(denom(jj,m))**(dfloat(m)/2.0d0)
 504  CONTINUE
      enddo
      !print *,'ENDSHERE',icc,idc,iec
      !print *,'sum',jjcount,sumdirect
      !print *,'Sum+Int=',sumdirect, jjcount

 507   format(i3,6f14.4)
 10    format(i10,1x,15G15.10)
 11    format(i5,i5,1x,15G15.10)

      RETURN
      END SUBROUTINE LAT_SUM

!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      SUBROUTINE LATSD(tvec,rr0,rcellw,ncellw,sumdirect)
      IMPLICIT REAL*8 (a-h,o-z)
      PARAMETER (LM=06)
      PARAMETER (MXPOISS=100)
      PARAMETER (LMX=3*LM)
      PARAMETER (MAX_ANG=((LMX+1)*(LMX+2))/2)
      PARAMETER (MPX=MAX_ANG)
      PARAMETER (MXPR=MXPOISS)

      DIMENSION SMULTI(100),ILM(100)
      DIMENSION tvec(3,3),rtvec(3),rr0(3),ccenter(3),rcellw(3,ncellw)
      DIMENSION denom(10000),unitang(3,10000),POLY(MPX,(LM+1)**2)

      CHARACTER(8) :: Date
      CHARACTER(10) :: Time1,Time2
      CHARACTER(5) :: Zone
      INTEGER :: Values(8)

      pi=4.0D0*atan(1.0D0)

      ccenter = 0.0d0
      rr0(1)=rr0(1)-ccenter(1)
      rr0(2)=rr0(2)-ccenter(2)
      rr0(3)=rr0(3)-ccenter(3)

      amaxc=0.0d0
      do ii=1,ncellw
         amaxc=max(amaxc,abs(rcellw(1,ii)))
      enddo
      print *,'amaxc',amaxc

      rtvec   = 0.0D0
      denom   = 0.0D0
      unitang = 0.0D0

      OPEN(299,FILE='SMULTI',FORM='formatted',STATUS='OLD')
      READ(299,*)ILMTOT
      DO kk=1,ILMTOT
      READ(299,*)ILM(kk),SMULTI(kk)
      ENDDO
      CLOSE(299)

!      do ii=1,ILMTOT
!         print *,'smulti00',ILM(ii),SMULTI(ii)
!      enddo

      iicount=0
      DO  n1=-30,30
         xnn=dfloat(n1)
! Removing cells considered in Wannier Domain 
!=================================================================
         if(abs(xnn).le.amaxc)goto 10
         iicount=iicount+1
         DO ix=1,3
         rtvec(ix)=rr0(ix)-xnn*TVEC(ix,1)
         ENDDO
         denom(iicount)=dsqrt(rtvec(1)**2+rtvec(2)**2+rtvec(3)**2)
         DO ix=1,3
         unitang(ix,iicount)=rtvec(ix)/(denom(iicount))
         ENDDO
  10  CONTINUE
      ENDDO

!      print *,'iicount',iicount

      POLY=0.0D0
      CALL HARMONICS(MPX,iicount,LM,unitang,POLY,NPOLY)

! NPOLY should be equal to ILMTOT
      sumdirect=0.0d0
      DO ii=1,iicount
      DO kk=1,NPOLY
         mm=ILM(kk)
         m=ceiling(dsqrt(dfloat(ILM(kk))))
         sumdirect=sumdirect+(SMULTI(kk)*POLY(ii,kk))
     &             /((denom(ii))**(dfloat(m)))
      ENDDO
      ENDDO

      RETURN
      END SUBROUTINE LATSD

!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      SUBROUTINE LAT_DIRECT(rr0,sumdirect)
      INCLUDE 'PARAMS'
      INCLUDE 'commons.inc'
       PARAMETER (LM=06)
       PARAMETER (LMX=3*LM)
       PARAMETER (MAX_ANG=((LMX+1)*(LMX+2))/2)
       PARAMETER (MPX=MAX_ANG)
       PARAMETER (MXPR=MXPOISS)
       PARAMETER (ILAT=20)
       PARAMETER (IGRID=(2*ILAT+1)**3)
       DIMENSION angrid(IGRID),ANGR(3,1000)
       DIMENSION rtvec(3),rr0(3),ccenter(3),tvecw(3,IGRID)
       DIMENSION denom(IGRID),unitang(3,IGRID),POLY(IGRID,(LM+1)**2)
       DIMENSION SMULTIR(100),POLYR(MPX,(LM+1)**2),vin(3),vout(3)
       COMMON/PTRANS/TVEC(3,3),atheta,iraxis,
     &               NXYZ,ISHELLX,ISHELLY,ISHELLZ       
       COMMON/COULDW/NANG,ANG(3,MAX_ANG),DOMEGANG(MAX_ANG),RAD_SPH
       COMMON/TMP1/COULOMB(MAX_PTS),RHOG(MAX_PTS,NVGRAD,MXSPN)
       COMMON/REPLICA/NMSHW,NDCELL,MMSH1,NATMS,
     &       ZWATM(1000),RATMS(3,1000),ITSCF
      COMMON/LMULTI/ILMTOTW,ILMW(100),SMULTIW(100)
      Character(8) :: Date
      Character(10) :: Time1,Time2
      Character(5) :: Zone
      Integer :: Values(8)
      LOGICAL FTVEC
      DATA FTVEC/.TRUE./

       pi=4.0D0*atan(1.0D0)

       ccenter = 0.0d0
       rr0(1)=rr0(1)-ccenter(1)
       rr0(2)=rr0(2)-ccenter(2)
       rr0(3)=rr0(3)-ccenter(3)
 
       rtvec   =0.0d0
       denom   =0.0d0
       unitang =0.0d0
       angrid = 0.0d0

C From negative to positive values of n1,n2, and n3

       if(FTVEC)then
       open(30,file="NGRID")
       rewind(30)
       read(30,*)iicount
       do ii=1,iicount
          read(30,*)n1,n2,n3,anglat
          angrid(ii)=anglat
          do ix=1,3
          tvecw(ix,ii)=(n1*TVEC(ix,1)+n2*TVEC(ix,2)+n3*TVEC(ix,3))
          enddo
       enddo
       FTVEC=.FALSE.
       close(30)
       endif

       do ii=1,iicount
          do ix=1,3
          rtvec(ix)=rr0(ix)-tvecw(ix,ii)
          enddo
          denom(ii)=dsqrt(rtvec(1)**2+rtvec(2)**2+rtvec(3)**2)
          do ix=1,3
          unitang(ix,ii)=rtvec(ix)/(denom(ii))
          enddo
       enddo

       POLY=0.0D0
       CALL HARMONICS(IGRID,iicount,LM,unitang,POLY,NPOLY)
    
!       OPEN(299,FILE='SMULTI',FORM='formatted',STATUS='OLD')
!       !READ(299,*)(ccenter(ii),ii=1,3)
!       READ(299,*)ILMTOT
!       DO kk=1,ILMTOT
!       READ(299,*)ILM(kk),SMULTI(kk)
!       ENDDO
!       CLOSE(299)

       sumdirect=0.0d0
       DO ii=1,iicount
          IF(abs(atheta).gt.0.1d-10)THEN
          DO ix=1,NANG
             vin(1)=ANG(1,ix)
             vin(2)=ANG(2,ix)
             vin(3)=ANG(3,ix)
             RANG  =angrid(ii)
             call rotvec(iraxis,RANG,vin,vout)
             ANGR(1,ix)=vout(1)
             ANGR(2,ix)=vout(2)
             ANGR(3,ix)=vout(3)
          ENDDO

          POLYR=0.0d0
          CALL HARMONICS(MPX,NANG,LM,ANGR,POLYR,NPOLYR)

          DO I=1,((LM+1)**2)
          SMULTIR(I)=0.0D0
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
             SMULTIR(I)=SMULTIR(I)+POLYR(kk,I)*DOMEGANG(kk)
     &                *(COULOMB(IPTS)+ZPOTTEMP)
             kk=kk+1
             ENDDO ! IPTS
             ll=ceiling(dsqrt(dfloat(I)))
             SMULTIR(I)=SMULTIR(I)*(rad_sph**dfloat(ll))
          ENDDO ! I

          ELSE
             DO ij=1,NPOLY
             SMULTIR(ij)=SMULTIW(ij)
             ENDDO
          ENDIF !ipax


          DO ik=1,NPOLY
             mm=ILMW(ik)
             m=ceiling(dsqrt(dfloat(ILMW(ik))))
             sumdirect=sumdirect+(SMULTIR(ik)*POLY(ii,ik))
     &                 /((denom(ii))**(dfloat(m)))
          ENDDO
         ENDDO

       RETURN
       END SUBROUTINE LAT_DIRECT

!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      ! The following is for periodic 2D system along X-direction.
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
       SUBROUTINE LAT_SUM2(rr0,rcellw,ncellw,sumdirect)
       IMPLICIT REAL*8 (a-h,o-z)
       PARAMETER (LM=06)
       PARAMETER (MXPOISS=100)
       PARAMETER (LMX=3*LM)
       PARAMETER (MAX_ANG=((LMX+1)*(LMX+2))/2)
       PARAMETER (MPX=MAX_ANG)
       PARAMETER (MXPR=MXPOISS)

       DIMENSION SMULTI(100),ILM(100)
       DIMENSION sump(10),alp(10),rtvec(3)
       DIMENSION rr0(3),ccenter(3)
       DIMENSION xn1(10),xn2(10),xn3(10),wn1(10),wn2(10),wn3(10)
       DIMENSION rcellw(3,10000),POLY(MPX,((LM+1)**2))

       ALLOCATABLE :: unitvec(:,:,:),unitvec2(:,:),denom(:,:),wtt(:,:)
       ALLOCATABLE :: xvec(:),yvec(:),zvec(:),wxvec(:),wyvec(:),wzvec(:)
       COMMON/PTRANS/TVEC(3,3),atheta,iraxis,
     &               NXYZ,ISHELLX,ISHELLY,ISHELLZ       

      Character(8) :: Date
      Character(10) :: Time1,Time2
      Character(5) :: Zone
      Integer :: Values(8)

      pi = 4.0D0*atan(1.0D0)

      ! Gauss quadrature nn-point formula
      nn = 5
      l=5

      do m=1,10
         sump(m)=0.0D0
         alp(m)=0.0D0
      enddo

      do n=1,l-1
        do m=2,10
          sump(m)=sump(m)+1.0D0/(dfloat(n)**dfloat(m))
        enddo
      enddo

      alp(2)=pi*pi/6.0D0-sump(2)
      alp(2)=1.0D0/alp(2)
      alp(3)=2.0D0*(1.20205690315959D0-sump(3))
      alp(3)=1.0D0/alp(3)**(1.0D0/2.0D0)
      alp(4)=3.0D0*(1.08232323371114D0-sump(4))
      alp(4)=1.0D0/alp(4)**(1.0D0/3.0D0)
! Zeta(5) is blowing and for now we approximate zeta(5) with zeta(4)
      alp(5)=alp(4)
!      alp(5)=4.0D0*(1.0369277551433699263-sump(5))
!      alp(5)=1.0D0/alp(5)**(1.0D0/4.0D0)
      alp(6)=5.0D0*(pi**6.0d0/945.0D0-sump(6))
      alp(6)=1.0D0/alp(6)**(1.0D0/5.0D0)

      do ii=7,10
        alp(ii)=alp(6)
      enddo

      igrid=(l**2+2*nn*l+2*nn**2)

      jgrid=4*igrid

      MM_MAX=(LM+1)**2

      ALLOCATE(wxvec(igrid))
      ALLOCATE(wyvec(igrid))
      ALLOCATE(wzvec(igrid))

      ALLOCATE(xvec(igrid))
      ALLOCATE(yvec(igrid))
      ALLOCATE(zvec(igrid))

      wxvec(:)=0.0d0
      wyvec(:)=0.0d0
      wzvec(:)=0.0d0

      xvec(:)=0.0d0
      yvec(:)=0.0d0
      zvec(:)=0.0d0

      do ik=1,10
        xn1(ik)=0.0d0
        xn2(ik)=0.0d0
        xn3(ik)=0.0d0
        wn1(ik)=0.0d0
        wn2(ik)=0.0d0
        wn3(ik)=0.0d0
      enddo

      OPEN(299,FILE='SMULTI',FORM='formatted',STATUS='OLD')
      REWIND(299)
      !READ(299,*)(ccenter(ii),ii=1,3)
      READ(299,*)ILMTOT

      do kk=1,ILMTOT
        READ(299,*)ILM(kk),SMULTI(kk)
      enddo
      CLOSE(299)

!      do kk=1,ILMTOT
!        print *,ILM(kk),SMULTI(kk)
!      enddo

      ccenter=0.0d0

      rr0(1)=rr0(1)-ccenter(1)
      rr0(2)=rr0(2)-ccenter(2)
      rr0(3)=rr0(3)-ccenter(3)

!      ishell=0
!      itotmultipole=7

!      do kk=1,ILMTOT
!Only considering first non-vanishing multipole moment
      !mm=ILM(1) 
      !m=ceiling(dsqrt(dfloat(mm)))
      !print *,'The m value=',m

      m=3 !hard coded
      iicount=0
      do n1=0,l-1
        do n2=0,l-1
           n3=0
          if(n1.eq.0.and.n2.eq.0.and.n3.eq.0)goto 800

          iicount=iicount+1

          wxvec(iicount)=1.0d0
          wyvec(iicount)=1.0d0
          wzvec(iicount)=1.0d0

          xvec(iicount)=dfloat(n1)
          yvec(iicount)=dfloat(n2)
          zvec(iicount)=dfloat(n3)

  800   continue
        enddo !n2 loop ends      
!==========================================================================

        dlimn20 = 0.0D0
        dlimn21 = 1.0D0/alp(m)

        call gaussp(dlimn20,dlimn21,nn,xn2,wn2)
        do j=1,nn
           iicount=iicount+1

           wxvec(iicount)=1.0d0
           wyvec(iicount)=wn2(j)/(xn2(j)**2.0d0)
           wzvec(iicount)=1.0d0

           xvec(iicount)=dfloat(n1)
           yvec(iicount)=1.0d0/xn2(j)
           zvec(iicount)=dfloat(n3)
        enddo
      enddo !n1 loop ends

      do n2=0,l-1
        limn10 = 0.0D0
        dlimn11 = 1.0D0/alp(m)

        call gaussp(dlimn10,dlimn11,nn,xn1,wn1)
        do i=1,nn
           iicount=iicount+1

           wxvec(iicount)=wn1(i)/(xn1(i)**2.0d0)
           wyvec(iicount)=1.0d0
           wzvec(iicount)=1.0d0

           xvec(iicount)=1.0d0/xn1(i)
           yvec(iicount)=dfloat(n2)
           zvec(iicount)=dfloat(n3)
        enddo
      enddo

!------------------------------------------------------------------------            
! Double Integration
!------------------------------------------------------------------------            

        dlimn10=0.0D0
        dlimn11=pi/4.0D0
        call gaussp(dlimn10,dlimn11,nn,xn1,wn1)
        do i=1,nn
           dlimn20=0.0D0
           dlimn21=1.0D0/(alp(m)*cos(xn1(i)))
           call gaussp(dlimn20,dlimn21,nn,xn2,wn2)
           do j=1,nn
              rr=xn2(j)
!========================================================================              
!             Term 1
              vv=rr*cos(xn1(i))
              ww=rr*sin(xn1(i))

              iicount=iicount+1

              wxvec(iicount)=wn1(i)*(dsqrt(rr)/vv**2.0d0)
              wyvec(iicount)=wn2(j)*(dsqrt(rr)/ww**2.0d0)
              wzvec(iicount)=1.0d0

              xvec(iicount)=1.0d0/vv
              yvec(iicount)=1.0d0/ww
              zvec(iicount)=dfloat(n3)

!========================================================================              
!             Term 2
              vv=rr*sin(xn1(i))
              ww=rr*cos(xn1(i))

              iicount=iicount+1

              wxvec(iicount)=wn1(i)*(dsqrt(rr)/vv**2.0d0)
              wyvec(iicount)=wn2(j)*(dsqrt(rr)/ww**2.0d0)
              wzvec(iicount)=1.0d0

              xvec(iicount)=1.0d0/vv
              yvec(iicount)=1.0d0/ww
              zvec(iicount)=dfloat(n3)
           enddo
        enddo

!  505  continue ! m loop started at the top ends here 

!        print *,iicount
!        print *,''
!        do ii=1,iicount
!        print 506,'He',xvec(ii),yvec(ii),zvec(ii),
!     &    wxvec(ii),wyvec(ii),wzvec(ii)
!        enddo
!  506   format(a3,6f14.4)
!        stop

      ALLOCATE(unitvec(3,jgrid,LM+1))
      ALLOCATE(denom(jgrid,LM+1))
      ALLOCATE(wtt(jgrid,LM+1))

      unitvec(:,:,:) = 0.0d0
      wtt(:,:) = 0.0d0
      denom(:,:) = 0.0d0
      rtvec(:) = 0.0d0

      m=3 !hard coding m
      jjcount=0
      do ii=1,iicount

        xnn=xvec(ii)
        ynn=yvec(ii)
        znn=zvec(ii)
        
!1------------- 1 0 0------------------------------
        if(xnn.ne.0.0d0.and.ynn.eq.0.0d0.and.znn.eq.0.0d0)then

          xnn2=xnn
          ynn2=ynn
          znn2=znn

          do jj=1,ncellw
             if(xnn2.eq.rcellw(1,jj).and.ynn2.eq.rcellw(2,jj)
     &       .and.znn2.eq.rcellw(3,jj)) goto 998
          end do 

          jjcount=jjcount+1
          wtt(jjcount,m)=wxvec(ii)*wyvec(ii)*wzvec(ii)
          do ix=1,3
          rtvec(ix)=rr0(ix)-xnn2*TVEC(ix,1)-ynn2
     &         *TVEC(ix,2)-znn2*TVEC(ix,3)
          enddo
          denom(jjcount,m)=dsqrt(rtvec(1)**2+rtvec(2)**2+rtvec(3)**2)
          do ix=1,3
          unitvec(ix,jjcount,m)=rtvec(ix)/(denom(jjcount,m))
          enddo
          !print 507,jjcount,xnn2,ynn2,znn2,absrtvec,wtt(jjcount,m)

!1------------- -1 0 0------------------------------
          jjcount=jjcount+1

          xnn2=-xnn
          ynn2=ynn
          znn2=znn

          do jj=1,ncellw
             if(xnn2.eq.rcellw(1,jj).and.ynn2.eq.rcellw(2,jj)
     &       .and.znn2.eq.rcellw(3,jj)) goto 998
          end do 

          jjcount=jjcount+1
          wtt(jjcount,m)=wxvec(ii)*wyvec(ii)*wzvec(ii)
          do ix=1,3
          rtvec(ix)=rr0(ix)-xnn2*TVEC(ix,1)-ynn2
     &         *TVEC(ix,2)-znn2*TVEC(ix,3)
          enddo
          denom(jjcount,m)=dsqrt(rtvec(1)**2+rtvec(2)**2+rtvec(3)**2)
          do ix=1,3
          unitvec(ix,jjcount,m)=rtvec(ix)/(denom(jjcount,m))
          enddo

        endif

!1------------- 0 1 0------------------------------
        if(xnn.eq.0.0d0.and.ynn.ne.0.0d0.and.znn.eq.0.0d0)then

          xnn2=xnn
          ynn2=ynn
          znn2=znn
          do jj=1,ncellw
             if(xnn2.eq.rcellw(1,jj).and.ynn2.eq.rcellw(2,jj)
     &       .and.znn2.eq.rcellw(3,jj)) goto 998
          end do
      
          jjcount=jjcount+1
          wtt(jjcount,m)=wxvec(ii)*wyvec(ii)*wzvec(ii)
          do ix=1,3
          rtvec(ix)=rr0(ix)-xnn2*TVEC(ix,1)-ynn2
     &         *TVEC(ix,2)-znn2*TVEC(ix,3)
          enddo
          denom(jjcount,m)=dsqrt(rtvec(1)**2+rtvec(2)**2+rtvec(3)**2)
          do ix=1,3
          unitvec(ix,jjcount,m)=rtvec(ix)/(denom(jjcount,m))
          enddo


!1------------- 0 -1 0------------------------------

          xnn2=xnn
          ynn2=-ynn
          znn2=znn
          do jj=1,ncellw
             if(xnn2.eq.rcellw(1,jj).and.ynn2.eq.rcellw(2,jj)
     &       .and.znn2.eq.rcellw(3,jj)) goto 998
          end do
      
          jjcount=jjcount+1
          wtt(jjcount,m)=wxvec(ii)*wyvec(ii)*wzvec(ii)
          do ix=1,3
          rtvec(ix)=rr0(ix)-xnn2*TVEC(ix,1)-ynn2
     &         *TVEC(ix,2)-znn2*TVEC(ix,3)
          enddo
          denom(jjcount,m)=dsqrt(rtvec(1)**2+rtvec(2)**2+rtvec(3)**2)
          do ix=1,3
          unitvec(ix,jjcount,m)=rtvec(ix)/(denom(jjcount,m))
          enddo

        endif

!1------------- 1 1 0------------------------------

        if(xnn.ne.0.0d0.and.ynn.ne.0.0d0.and.znn.eq.0.0d0)then

          xnn2=xnn
          ynn2=ynn
          znn2=znn
          do jj=1,ncellw
             if(xnn2.eq.rcellw(1,jj).and.ynn2.eq.rcellw(2,jj)
     &       .and.znn2.eq.rcellw(3,jj)) goto 998
          end do
      
          jjcount=jjcount+1
          wtt(jjcount,m)=wxvec(ii)*wyvec(ii)*wzvec(ii)
          do ix=1,3
          rtvec(ix)=rr0(ix)-xnn2*TVEC(ix,1)-ynn2
     &         *TVEC(ix,2)-znn2*TVEC(ix,3)
          enddo
          denom(jjcount,m)=dsqrt(rtvec(1)**2+rtvec(2)**2+rtvec(3)**2)
          do ix=1,3
          unitvec(ix,jjcount,m)=rtvec(ix)/(denom(jjcount,m))
          enddo


!1------------- -1 1 0------------------------------

          xnn2=-xnn
          ynn2=ynn
          znn2=znn
          do jj=1,ncellw
             if(xnn2.eq.rcellw(1,jj).and.ynn2.eq.rcellw(2,jj)
     &       .and.znn2.eq.rcellw(3,jj)) goto 998
          end do
      
          jjcount=jjcount+1
          wtt(jjcount,m)=wxvec(ii)*wyvec(ii)*wzvec(ii)
          do ix=1,3
          rtvec(ix)=rr0(ix)-xnn2*TVEC(ix,1)-ynn2
     &         *TVEC(ix,2)-znn2*TVEC(ix,3)
          enddo
          denom(jjcount,m)=dsqrt(rtvec(1)**2+rtvec(2)**2+rtvec(3)**2)
          do ix=1,3
          unitvec(ix,jjcount,m)=rtvec(ix)/(denom(jjcount,m))
          enddo

!1------------- 1 -1 0------------------------------

          xnn2=xnn
          ynn2=-ynn
          znn2=znn
          do jj=1,ncellw
             if(xnn2.eq.rcellw(1,jj).and.ynn2.eq.rcellw(2,jj)
     &       .and.znn2.eq.rcellw(3,jj)) goto 998
          end do
      
          jjcount=jjcount+1
          wtt(jjcount,m)=wxvec(ii)*wyvec(ii)*wzvec(ii)
          do ix=1,3
          rtvec(ix)=rr0(ix)-xnn2*TVEC(ix,1)-ynn2
     &         *TVEC(ix,2)-znn2*TVEC(ix,3)
          enddo
          denom(jjcount,m)=dsqrt(rtvec(1)**2+rtvec(2)**2+rtvec(3)**2)
          do ix=1,3
          unitvec(ix,jjcount,m)=rtvec(ix)/(denom(jjcount,m))
          enddo

!1------------- -1 -1 0------------------------------

          xnn2=-xnn
          ynn2=-ynn
          znn2=znn
          do jj=1,ncellw
             if(xnn2.eq.rcellw(1,jj).and.ynn2.eq.rcellw(2,jj)
     &       .and.znn2.eq.rcellw(3,jj)) goto 998
          end do
      
          jjcount=jjcount+1
          wtt(jjcount,m)=wxvec(ii)*wyvec(ii)*wzvec(ii)
          do ix=1,3
          rtvec(ix)=rr0(ix)-xnn2*TVEC(ix,1)-ynn2
     &         *TVEC(ix,2)-znn2*TVEC(ix,3)
          enddo
          denom(jjcount,m)=dsqrt(rtvec(1)**2+rtvec(2)**2+rtvec(3)**2)
          do ix=1,3
          unitvec(ix,jjcount,m)=rtvec(ix)/(denom(jjcount,m))
          enddo

        endif

  998   continue
        enddo

        DEALLOCATE(xvec,yvec,zvec,wxvec,wyvec,wzvec)
! REARRAGE jjcount in increasing order of radius


!      do jj=1,jjcount
!        print *,(unitvec(ix,jj,m),ix=1,3)
!      enddo

!      if(igrid.eq.iicount)then
!        do m=2,itotmultipole

!        print *,'Total points i{n1 n2 n3}=',jjcount
        ALLOCATE(unitvec2(3,jjcount))

        do jj=1,jjcount
        unitvec2(1,jj)=unitvec(1,jj,m)
        unitvec2(2,jj)=unitvec(2,jj,m)
        unitvec2(3,jj)=unitvec(3,jj,m)
        enddo

!        do jj=1,jjcount
!          print *,(unitvec2(ix,jj),ix=1,3)
!        enddo

        CALL HARMONICS(MPX,jjcount,LM,unitvec2,POLY,NPOLY)

! Checking the output YLM generated above
!      print *,NPOLY
!      do i=1,NPOLY
!        do j=1,iicount
!          print *,POLY(j,i)
!        enddo
!      enddo

      sumdirect=0.0d0

      do jj=1,jjcount
         do kk=1,ILMTOT
           mm=ILM(kk)
           m=ceiling(dsqrt(dfloat(mm)))
           sumdirect=sumdirect+(SMULTI(kk)*POLY(jj,kk)*wtt(jj,3))
     &             /((denom(jj,3))**(dfloat(m)))
         enddo
      enddo
      !print *,'sum',jjcount,sumdirect
      !print *,'Sum+Int=',sumdirect, jjcount

      DEALLOCATE(unitvec,unitvec2,denom,wtt)

 507   format(i3,6f14.4)
 10    format(i10,1x,15G15.10)
 11    format(i5,i5,1x,15G15.10)

      RETURN
      END SUBROUTINE LAT_SUM2

!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
       SUBROUTINE LAT_SUM1(rr0,sumdirect)
       IMPLICIT REAL*8 (a-h,o-z)
       PARAMETER (LM=06)
       PARAMETER (MXPOISS=100)
       PARAMETER (LMX=3*LM)
       PARAMETER (MAX_ANG=((LMX+1)*(LMX+2))/2)
       PARAMETER (MPX=MAX_ANG)
       PARAMETER (MXPR=MXPOISS)

       DIMENSION SMULTI(100),ILM(100)
       DIMENSION sump(10),alp(10),rtvec(3)
       DIMENSION rr0(3),ccenter(3)
       DIMENSION xn1(10),xn2(10),xn3(10),wn1(10),wn2(10),wn3(10)
       DIMENSION rcellw(3,10000),POLY(MPX,(LM+1)**2)

       ALLOCATABLE :: unitvec(:,:),unitvec2(:,:),denom(:),wtt(:)
       ALLOCATABLE :: xvec(:),yvec(:),zvec(:),wxvec(:),wyvec(:),wzvec(:)
       COMMON/PTRANS/TVEC(3,3),atheta,iraxis,
     &               NXYZ,ISHELLX,ISHELLY,ISHELLZ       

       Character(8) :: Date
       Character(10) :: Time1,Time2
       Character(5) :: Zone
       Integer :: Values(8)
 
       pi = 4.0D0*atan(1.0D0)
 
       !print *,'Total cell in wannier domain',ncellw
       ! Gauss quadrature nn-point formula
       nn = 5
       l=7
 
       do m=1,10
         sump(m)=0.0D0
         alp(m)=0.0D0
       enddo
 
       do n=1,l-1
       do m=2,10
         sump(m)=sump(m)+1.0D0/(dfloat(n)**dfloat(m))
       enddo
       enddo

!===========================================================================
      alp(2)=pi*pi/6.0D0-sump(2)
      alp(2)=1.0D0/alp(2)
      alp(3)=2.0D0*(1.20205690315959D0-sump(3))
      alp(3)=1.0D0/alp(3)**(1.0D0/2.0D0)
      alp(4)=3.0D0*(1.08232323371114D0-sump(4))
      alp(4)=1.0D0/alp(4)**(1.0D0/3.0D0)
! Zeta(5) is blowing and for now we approximate zeta(5) with zeta(4)
      alp(5)=alp(4)
!      alp(5)=4.0D0*(1.0369277551433699263-sump(5))
!      alp(5)=1.0D0/alp(5)**(1.0D0/4.0D0)
      alp(6)=5.0D0*(pi**6.0d0/945.0D0-sump(6))
      alp(6)=1.0D0/alp(6)**(1.0D0/5.0D0)

      do ii=7,10
        alp(ii)=alp(6)
      enddo
!===========================================================================

      igrid=2*(l+nn)

      MM_MAX=(LM+1)**2

      ALLOCATE(wxvec(igrid))
      ALLOCATE(wyvec(igrid))
      ALLOCATE(wzvec(igrid))

      ALLOCATE(xvec(igrid))
      ALLOCATE(yvec(igrid))
      ALLOCATE(zvec(igrid))

      wxvec(:)=0.0d0
      wyvec(:)=0.0d0
      wzvec(:)=0.0d0

      xvec(:)=0.0d0
      yvec(:)=0.0d0
      zvec(:)=0.0d0

      do ik=1,10
        xn1(ik)=0.0d0
        xn2(ik)=0.0d0
        xn3(ik)=0.0d0
        wn1(ik)=0.0d0
        wn2(ik)=0.0d0
        wn3(ik)=0.0d0
      enddo

      OPEN(299,FILE='SMULTI',FORM='formatted',STATUS='OLD')
      READ(299,*)ILMTOT
      DO kk=1,ILMTOT
      READ(299,*)ILM(kk),SMULTI(kk)
      ENDDO
      CLOSE(299)

      ccenter=0.0d0

      rr0(1)=rr0(1)-ccenter(1)
      rr0(2)=rr0(2)-ccenter(2)
      rr0(3)=rr0(3)-ccenter(3)

! Getting grid points
!=================================================
      iicount=0
      do n1=0,l-1
         n3=0
         n2=0
        if(n1.eq.0.and.n2.eq.0.and.n3.eq.0)goto 800
        iicount=iicount+1
        wxvec(iicount)=1.0d0
        wyvec(iicount)=1.0d0
        wzvec(iicount)=1.0d0
        xvec(iicount)=dfloat(n1)
        yvec(iicount)=dfloat(n2)
        zvec(iicount)=dfloat(n3)
 800  continue
      enddo

      ! hard coded for m
      m=3
!      do m=1,6
       limn10 = 0.0D0
       dlimn11 = 1.0D0/alp(m)

       call gaussp(dlimn10,dlimn11,nn,xn1,wn1)

      do i=1,nn
         iicount=iicount+1

         wxvec(iicount)=wn1(i)/(xn1(i)**2.0d0)
         wyvec(iicount)=1.0d0
         wzvec(iicount)=1.0d0

         xvec(iicount)=1.0d0/xn1(i)
         yvec(iicount)=dfloat(n2)
         zvec(iicount)=dfloat(n3)
      enddo
!      enddo

!=================================================
!      print *,'xvec(ii),yvec(ii),zvec(ii),wxvec(ii),wyvec(ii),wzvec(ii)'
!      do ii=1,iicount
!       print *,xvec(ii),yvec(ii),zvec(ii),wxvec(ii),wyvec(ii),wzvec(ii)
!      enddo
!      call stopit

      ALLOCATE(unitvec(3,igrid))
      ALLOCATE(denom(igrid))
      ALLOCATE(wtt(igrid))

      unitvec(:,:) = 0.0d0
      wtt(:) = 0.0d0
      denom(:) = 0.0d0
      rtvec(:) = 0.0d0

      jjcount=0

      do ii=1,iicount
        xnn=xvec(ii)
        ynn=yvec(ii)
        znn=zvec(ii)

!1------------- 1 0 0------------------------------
        if(xnn.ne.0.0d0.and.ynn.eq.0.0d0.and.znn.eq.0.0d0)then

          xnn2=xnn
          ynn2=ynn
          znn2=znn

        do jj=1,ncellw
          if(xnn2.eq.rcellw(1,jj).and.ynn2.eq.rcellw(2,jj).and.
     &       znn2.eq.rcellw(3,jj))goto 998
        end do 

          jjcount=jjcount+1
          wtt(jjcount)=wxvec(ii)*wyvec(ii)*wzvec(ii)
          do ix=1,3
          rtvec(ix)=rr0(ix)-xnn2*TVEC(ix,1)-ynn2
     &         *TVEC(ix,2)-znn2*TVEC(ix,3)
          enddo
          denom(jjcount)=dsqrt(rtvec(1)**2+rtvec(2)**2+rtvec(3)**2)
          do ix=1,3
          unitvec(ix,jjcount)=rtvec(ix)/(denom(jjcount))
          enddo
          !print 507,jjcount,xnn2,ynn2,znn2,absrtvec,wtt(jjcount,m)

!1------------- -1 0 0------------------------------

          xnn2=-xnn
          ynn2=ynn
          znn2=znn

        do jj=1,ncellw
          if(xnn2.eq.rcellw(1,jj).and.ynn2.eq.rcellw(2,jj).and.
     &       znn2.eq.rcellw(3,jj))goto 998
        end do 
          jjcount=jjcount+1
          wtt(jjcount)=wxvec(ii)*wyvec(ii)*wzvec(ii)
          do ix=1,3
          rtvec(ix)=rr0(ix)-xnn2*TVEC(ix,1)-ynn2
     &         *TVEC(ix,2)-znn2*TVEC(ix,3)
          enddo
          denom(jjcount)=dsqrt(rtvec(1)**2+rtvec(2)**2+rtvec(3)**2)
          do ix=1,3
          unitvec(ix,jjcount)=rtvec(ix)/(denom(jjcount))
          enddo
        endif
  998   continue
        enddo
        DEALLOCATE(xvec,yvec,zvec,wxvec,wyvec,wzvec)

! REARRAGE jjcount in increasing order of radius


!      do jj=1,jjcount
!        print *,(unitvec(ix,jj,m),ix=1,3)
!      enddo

!      if(igrid.eq.iicount)then
!        do m=2,itotmultipole

!        print *,'Total points i{n1 n2 n3}=',jjcount
!        ALLOCATE(unitvec2(3,jjcount))
!
!        do jj=1,jjcount
!        unitvec2(1,jj)=unitvec(1,jj,m)
!        unitvec2(2,jj)=unitvec(2,jj,m)
!        unitvec2(3,jj)=unitvec(3,jj,m)
!        enddo
!
!        do jj=1,jjcount
!          print *,(unitvec2(ix,jj),ix=1,3)
!        enddo

        CALL HARMONICS(MPX,jjcount,LM,unitvec,POLY,NPOLY)

! Checking the output YLM generated above
!      print *,NPOLY
!      do i=1,NPOLY
!        do j=1,iicount
!          print *,POLY(j,i)
!        enddo
!      enddo

      sumdirect=0.0d0
      do jj=1,jjcount
         do kk=1,NPOLY
           mm=ILM(kk)
           m=ceiling(dsqrt(dfloat(mm)))
           sumdirect=sumdirect+(SMULTI(kk)*POLY(jj,kk)*wtt(jj))
     &             /(denom(jj))**(dfloat(m))
         enddo
      enddo
      !print *,'sum',jjcount,sumdirect
      !print *,'Sum+Int=',sumdirect, jjcount

      DEALLOCATE(unitvec,denom,wtt)

! 507   format(i3,6f14.4)
! 10    format(i10,1x,15G15.10)
! 11    format(i5,i5,1x,15G15.10)

      RETURN
      END SUBROUTINE LAT_SUM1

!!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!      ! Check SMULTI
!!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!!      SUBROUTINE TEST_SMULTI(TVEC,r0,atom2)
!!      IMPLICIT REAL*8 (a-h,o-z)
!!      DIMENSION TVEC(3,3),r0(3),atom2(8,3)
!!      DIMENSION vec(3),atomvec(3),vect_nn(3)
!!      DIMENSION unitvec(3,1),POLY(1,49)
!!      DIMENSION ILM(49),SMULTI(49)
!!
!!      PARAMETER (LM=06)
!!
!!      coupot=0.0d0
!!
!!
!!      do iat=1,8
!!         atomvec(:)=atom2(iat,:)
!!         vec(:)=r0(:)-atomvec(:)
!!         denom1=dsqrt(sum(vec(:)*vec(:)))
!!        if(iat.le.4)then
!!          coupot=coupot+1.0d0/denom1
!!        else
!!          coupot=coupot-1.0d0/denom1
!!        endif
!!       enddo
!!
!!      print *,'Potential Direct=',coupot
!!      
!!
!!      dist_vec=dsqrt(sum(r0(:)*r0(:)))
!!
!!      do ix=1,3
!!      unitvec(ix,1)=r0(ix)/dist_vec
!!      enddo
!!
!!      CALL HARMONICS(1,1,LM,unitvec,POLY,NPOLY)
!!
!!!      do jj=1,49
!!!      print *,POLY(1,jj)
!!!      enddo
!!
!!      OPEN(299,FILE='SMULTI',FORM='formatted',STATUS='OLD')
!!      REWIND(299)
!!      READ(299,*)ILMTOT
!!      do kk=1,ILMTOT
!!         READ(299,*) ILM(kk),SMULTI(kk)
!!      enddo
!!
!!      sumdirect=0.0d0
!!
!!         do mm=1,ILMTOT
!!           lms=ILM(mm) 
!!           m=ceiling(dsqrt(dfloat(lms)))
!!           sumdirect=sumdirect+(SMULTI(mm)*POLY(1,lms))/
!!     &            (dist_vec)**(dfloat(m))
!!         enddo
!!
!!      print *,'Potential SMULTI=',sumdirect
!!      CLOSE(299)
!!      RETURN
!!
!!      END SUBROUTINE
!
!c
!c  routine to decompose the coulomb potential of the central cell into multipole moments (spherical harmonics)
!c
