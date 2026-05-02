      PROGRAM PERIODIC_3D
      IMPLICIT REAL*8 (a-h,o-z)
      DIMENSION ttmat(3, 3),rr0(3)
      ALLOCATABLE :: atoms(:,:)

      pi = 4.0D0*atan(1.0D0)

      natoms=2

      open(99,file='lattice.data', form='formatted', status='old')
      READ(99,*)
      READ(99,*)
      do i=1,3
        READ(99,*)(ttmat(i,j),j=1,3)
      enddo
      
      READ(99,*)
      READ(99,*)
      READ(99,*)

      ALLOCATE(atoms(natoms,3))

      do iatom=1,natoms
          read(99,*) (atoms(iatom,j), j = 1, 3)
      enddo

      rr0(1) = 0.1D0
      rr0(2) = 0.1D0
      rr0(3) = 0.1D0

!      CALL TEST_SMULTI(ttmat,rr0,atom1)
      CALL LS_DIRECT(ttmat,rr0,atoms,natoms)
      CALL LAT_SUM(ttmat,rr0,atoms,natoms)

      DEALLOCATE(atoms)

      END PROGRAM PERIODIC_3D

!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      ! Check SMULTI
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      SUBROUTINE TEST_SMULTI(tmat,r0,atom2)
      IMPLICIT REAL*8 (a-h,o-z)
      DIMENSION tmat(3,3),r0(3),atom2(8,3)
      DIMENSION vec(3),atomvec(3),vect_nn(3)
      DIMENSION unitvec(3,1),POLY(1,49)
      DIMENSION ILM(49),SMULTI(49)

      PARAMETER (LM=06)

      coupot=0.0d0


      do iat=1,8
         atomvec(:)=atom2(iat,:)
         vec(:)=r0(:)-atomvec(:)
         denom1=dsqrt(sum(vec(:)*vec(:)))
        if(iat.le.4)then
          coupot=coupot+1.0d0/denom1
        else
          coupot=coupot-1.0d0/denom1
        endif
       enddo

      print *,'Potential Direct=',coupot
      

      dist_vec=dsqrt(sum(r0(:)*r0(:)))

      do ix=1,3
      unitvec(ix,1)=r0(ix)/dist_vec
      enddo

      CALL HARMONICS(1,1,LM,unitvec,POLY,NPOLY)

!      do jj=1,49
!      print *,POLY(1,jj)
!      enddo

      OPEN(299,FILE='SMULTI',FORM='formatted',STATUS='OLD')
      REWIND(299)
      READ(299,*)ILMTOT
      do kk=1,ILMTOT
         READ(299,*) ILM(kk),SMULTI(kk)
      enddo

      sumdirect=0.0d0

         do mm=1,ILMTOT
           lms=ILM(mm) 
           m=ceiling(dsqrt(dfloat(lms)))
           sumdirect=sumdirect+(SMULTI(mm)*POLY(1,lms))/
     &            (dist_vec)**(dfloat(m))
         enddo

      print *,'Potential SMULTI=',sumdirect
      CLOSE(299)

      END SUBROUTINE
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      SUBROUTINE LS_DIRECT(tmat,r0,atomloc,natom)
      IMPLICIT REAL*8 (a-h,o-z)
      DIMENSION tmat(3, 3),r0(3)

      PARAMETER (LM=08)
      PARAMETER (MXPOISS=100)
      PARAMETER (LMX=3*LM)
      PARAMETER (MAX_ANG=((LMX+1)*(LMX+2))/2)
      PARAMETER (MPX=MAX_ANG)
      PARAMETER (MXPR=MXPOISS)
      DIMENSION SMULTI((LM+1)**2)

      DIMENSION rtvec(3),r0d(3),ccenter(3),atomloc(natom,3)
      DIMENSION vect_nn(3),atomvec(3),vec(3)
      DIMENSION hisym(3,10),hincrement(10,3)

      ALLOCATABLE :: denom(:),unitang(:,:),POLY(:,:),ILM(:)

      Character(8) :: Date
      Character(10) :: Time1,Time2
      Character(5) :: Zone
      Integer :: Values(8)

!      print *,r0d(1),r0d(2),r0d(3)

      pi=4.0D0*atan(1.0D0)


!      open(399,file='output.dat',status='unknown',form='formatted')
!      rewind(399)
!
!      open(101,file='HISYMM',form='formatted',status='unknown')
!      rewind(101)
!      read(101,*)ihigrid,nhisymm
!      do ii=1,nhisymm ! Total number of high symm points
!        read(101,*)(hisym(ix,ii),ix=1,3)
!      enddo
!
!      do ii=1,nhisymm
!        do ix=1,3
!         hincrement(ix,ii)=(hisym(ix,ii+1)-hisym(ix,ii))/dfloat(ihigrid)
!        enddo
!      enddo
!
!      rmag_ad=0.0d0

      ilat=150
      igrid=(2*ilat+1)**3
 
      ALLOCATE(denom(igrid))
      ALLOCATE(unitang(3,igrid))

!      r0d(1)=0.0d0
!      r0d(2)=0.0d0
!      r0d(3)=0.0d0

!      do 801 inn=1,nhisymm-1
!        r0(:)=hisym(:,inn)

       OPEN(299,FILE='SMULTI',FORM='formatted',STATUS='OLD')
       REWIND(299)
       READ(299,*)(ccenter(ii),ii=1,3)
       READ(299,*)ILMTOT

       ALLOCATE(ILM(ILMTOT))

       do kk=1,ILMTOT
          READ(299,*)ILM(kk),SMULTI(kk)
       enddo

       r0d(1)=r0(1)-ccenter(1)
       r0d(2)=r0(2)-ccenter(2)
       r0d(3)=r0(3)-ccenter(3)

      print *,r0d
!       do 802 isym=1,ihigrid

       rtvec(:)=0.0d0
       denom(:)=0.0d0
       unitang(3,:)=0.0d0

       iicount=0
       ishell=3
C From negative to positive values of n1,n2, and n3
      do n1=-ilat,ilat
         do n2=-ilat,ilat
            do n3=-ilat,ilat

            if(n1.eq.0.and.n2.eq.0.and.n3.eq.0)goto 800
            if(abs(n1).le.ishell.and.abs(n2).le.ishell
     &           .and.abs(n3).le.ishell)goto 800
              n11=dfloat(n1)
              n22=dfloat(n2)
              n33=dfloat(n3)
              iicount=iicount+1
              rtvec(:)=(r0d(:)-n11*tmat(1,:)-n22
     &              *tmat(2,:)-n33*tmat(3,:))
              absrtvec=dsqrt(sum(rtvec(:)*rtvec(:)))
            do ix=1,3
              unitang(ix,iicount)=rtvec(ix)/absrtvec
            enddo
              denom(iicount)=(sum(rtvec(:)*rtvec(:)))
            !print *,iicount,n11,n22,n33,absrtvec
  800       continue
            enddo
         enddo
      enddo
      
      ALLOCATE(POLY(iicount,(LM+1)**2))
      POLY=0.0D0

      CALL HARMONICS(iicount,iicount,LM,unitang,POLY,NPOLY)

      sumdirect=0.0d0

      do ii=1,iicount
         do mm=1,ILMTOT
           lms=ILM(mm) 
           m=ceiling(dsqrt(dfloat(lms)))
           sumdirect=sumdirect+(SMULTI(mm)*POLY(ii,lms))/
     &            ((denom(ii))**(dfloat(m)/2.0d0))
         enddo
      enddo

!      do iat=1,natom
!         atomvec(:)=atom2(iat,:)
!         vec(:)=r0(:)-atomvec(:)
!         denom1=dsqrt(sum(vec(:)*vec(:)))
!        if(iat.le.natom/2)then
!          sumdirect=sumdirect+1.0d0/denom1
!        else
!          sumdirect=sumdirect-1.0d0/denom1
!        endif
!      enddo

!      rmag=dsqrt(abs(sum(r0(:)*r0(:))))
      !print *,rmag,rmag_ad,sumdirect
      print *,'Direct sum=',sumdirect

!      rold1=r0(1)
!      rold2=r0(2)
!      rold3=r0(3)
!
!      r0(1)=r0(1)+hincrement(1,inn)
!      r0(2)=r0(2)+hincrement(2,inn)
!      r0(3)=r0(3)+hincrement(3,inn)
!
!      addr=dsqrt((r0(1)-rold1)**2.0d0+(r0(2)-rold2)**2.0d0
!     &    +(r0(3)-rold3)**2.0d0)
!      rmag_ad=rmag_ad+addr

      DEALLOCATE(POLY)

!  802 continue
      DEALLOCATE(ILM)
!      CLOSE(299)
!  801 continue

      DEALLOCATE(denom,unitang)

 508   format(i3,4f14.4)
!      CLOSE(399)
!      CLOSE(101)

      END SUBROUTINE LS_DIRECT

!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      SUBROUTINE LAT_SUM(tmat,r0,atomloc,natom)
       IMPLICIT REAL*8 (a-h,o-z)
       PARAMETER (LM=08)
       PARAMETER (MXPOISS=100)
       PARAMETER (LMX=3*LM)
       PARAMETER (MAX_ANG=((LMX+1)*(LMX+2))/2)
       PARAMETER (MPX=MAX_ANG)
       PARAMETER (MXPR=MXPOISS)

       DIMENSION SMULTI((LM+1)**2)
       DIMENSION sump(10),alp(10),rtvec(3)
       DIMENSION tmat(3,3),r0(3),atomloc(natom,3),ccenter(3)
       DIMENSION xn1(10),xn2(10),xn3(10),wn1(10),wn2(10),wn3(10)

       ALLOCATABLE :: unitvec(:,:,:),unitvec2(:,:),denom(:,:),wtt(:,:)
       ALLOCATABLE :: xvec(:),yvec(:),zvec(:),wxvec(:),wyvec(:),wzvec(:)
       ALLOCATABLE :: POLY(:,:)
       ALLOCATABLE :: ILM(:)

      Character(8) :: Date
      Character(10) :: Time1,Time2
      Character(5) :: Zone
      Integer :: Values(8)

      pi = 4.0D0*atan(1.0D0)

      !print *,'Lattice sum begins!'
      ! Gauss quadrature nn-point formula
      nn = 5
      l=8

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

!      do ii=1,10
!      print *,ii,alp(ii)
!      enddo

      igrid=(l**3+3*nn*(l**2)+
     &      6*(nn**2)*l+6*nn**3)

      jgrid=8*igrid

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
      READ(299,*)(ccenter(ii),ii=1,3)
      READ(299,*)ILMTOT
      ALLOCATE(ILM(ILMTOT))

      do kk=1,ILMTOT
         READ(299,*) ILM(kk),SMULTI(kk)
      enddo

!      do kk=1,ILMTOT
!         print *,ILM(kk), SMULTI(kk)
!      enddo
!      stop

      r0(1)=r0(1)-ccenter(1)
      r0(2)=r0(2)-ccenter(2)
      r0(3)=r0(3)-ccenter(3)
      
      print *,r0

      ishell=3
      itotmultipole=7

!      do kk=1,ILMTOT
      mm=ILM(2)
      m=ceiling(dsqrt(dfloat(mm)))
      iicount=0
      !  print *,'LM, l+1 value',mm,m
!        do m=2,itotmultipole
        do n1=0,l-1
          do n2=0,l-1
            do n3=0,l-1

              if(n1.eq.0.and.n2.eq.0.and.n3.eq.0)goto 800
              if(n1.le.ishell.and.n2.le.ishell.and.n3
     &          .le.ishell)goto 800

               iicount=iicount+1

               wxvec(iicount)=1.0d0
               wyvec(iicount)=1.0d0
               wzvec(iicount)=1.0d0
               xvec(iicount)=dfloat(n1)
               yvec(iicount)=dfloat(n2)
               zvec(iicount)=dfloat(n3)

  800        continue
            enddo ! n3 loop ends
!==========================================================================
            !goto 504

            dlimn30 = 0.0D0
            dlimn31 = 1.0D0/alp(m)

            call gaussp(dlimn30,dlimn31,nn,xn3,wn3)
            do k=1,nn

               iicount=iicount+1

               wxvec(iicount)=1.0d0
               wyvec(iicount)=1.0d0
               wzvec(iicount)=wn3(k)/(xn3(k)**2.0d0)

               xvec(iicount)=dfloat(n1)
               yvec(iicount)=dfloat(n2)
               zvec(iicount)=1.0d0/xn3(k)
            enddo

  504   continue
        enddo !n2 loop ends      
!------------------------------------------------------------------------            
        !goto 503

        dlimn20=0.0D0
        dlimn21=pi/4.0D0
        call gaussp(dlimn20,dlimn21,nn,xn2,wn2)
        do j=1,nn
           dlimn30=0.0D0
           dlimn31=1.0D0/(alp(m)*cos(xn2(j)))
           call gaussp(dlimn30,dlimn31,nn,xn3,wn3)
           do k=1,nn

              rr=xn3(k)
              vv=rr*cos(xn2(j))
              ww=rr*sin(xn2(j))

              iicount=iicount+1

              wxvec(iicount)=1.0d0
              wyvec(iicount)=wn2(j)*(dsqrt(rr)/vv**2.0d0)
              wzvec(iicount)=wn3(k)*(dsqrt(rr)/ww**2.0d0)
              xvec(iicount)=dfloat(n1)
              yvec(iicount)=1.0d0/vv
              zvec(iicount)=1.0d0/ww

!========================================================================              
!             Term 2
              vv=rr*sin(xn2(j))
              ww=rr*cos(xn2(j))

              iicount=iicount+1

              wxvec(iicount)=1.0d0
              wyvec(iicount)=wn2(j)*(dsqrt(rr)/vv**2.0d0)
              wzvec(iicount)=wn3(k)*(dsqrt(rr)/ww**2.0d0)
              xvec(iicount)=dfloat(n1)
              yvec(iicount)=1.0d0/vv
              zvec(iicount)=1.0d0/ww
           enddo
        enddo
 503  continue
      enddo ! n1 loop ends.       
!=====================================================================              
      !goto 505 !goes to the enddo of m

      do n3=0,l-1
        do n1=0,l-1
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
         enddo ! n1 loop ends
!=======================================================================

         dlimn20 = 0.0D0
         dlimn21 = pi/4.0D0
         call gaussp(dlimn20,dlimn21,nn,xn2,wn2)
         do j=1,nn
            dlimn10 = 0.0D0
            dlimn11 = 1.0D0/alp(m)*cos(xn2(j))
            call gaussp(dlimn10,dlimn11,nn,xn1,wn1)
            do i=1,nn
              rr=xn1(i)
              uu=rr*cos(xn2(j))
              vv=rr*sin(xn2(j))
              iicount=iicount+1

              wxvec(iicount)=wn1(i)*(dsqrt(rr)/uu**2.0d0)
              wyvec(iicount)=wn2(j)*(dsqrt(rr)/vv**2.0d0)
              wzvec(iicount)=1.0d0
              xvec(iicount)=1.0d0/uu
              yvec(iicount)=1.0d0/vv
              zvec(iicount)=dfloat(n3)

!=====================================================================
              ! PART 2

              uu=rr*sin(xn2(j))
              vv=rr*cos(xn2(j))

              iicount=iicount+1

              wxvec(iicount)=wn1(i)*(dsqrt(rr)/uu**2.0d0)
              wyvec(iicount)=wn2(j)*(dsqrt(rr)/vv**2.0d0)
              wzvec(iicount)=1.0d0
              xvec(iicount)=1.0d0/uu
              yvec(iicount)=1.0d0/vv
              zvec(iicount)=dfloat(n3)

            enddo
         enddo
        enddo ! n3 loop      
!========================================================================

        do n2=0,l-1
           do n3=0,l-1
              dlimn10 = 0.0D0
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
           enddo ! n3 loop ends

!====================================================================================

           dlimn10 = 0.0D0
           dlimn11 = pi/4.0D0
           call gaussp(dlimn10,dlimn11,nn,xn1,wn1)
           do i=1,nn
              dlimn30 = 0.0D0
              dlimn31 = 1.0D0/(alp(m)*cos(xn1(i)))
              call gaussp(dlimn30,dlimn31,nn,xn3,wn3)
              do k=1,nn
                 rr=xn3(k)
                 uu=rr*cos(xn1(i))
                 ww=rr*sin(xn1(i))

                 iicount=iicount+1

                 wxvec(iicount)=wn1(i)*(dsqrt(rr)/uu**2.0d0)
                 wyvec(iicount)=1.0d0
                 wzvec(iicount)=wn3(k)*(dsqrt(rr)/ww**2.0d0)
                 xvec(iicount)=1.0d0/uu
                 yvec(iicount)=dfloat(n2)
                 zvec(iicount)=1.0d0/ww

!=========================================================================
                 uu=rr*sin(xn1(i))
                 ww=rr*cos(xn1(i))

                 iicount=iicount+1

                 wxvec(iicount)=wn1(i)*(dsqrt(rr)/uu**2.0d0)
                 wyvec(iicount)=1.0d0
                 wzvec(iicount)=wn3(k)*(dsqrt(rr)/ww**2.0d0)
                 xvec(iicount)=1.0d0/uu
                 yvec(iicount)=dfloat(n2)
                 zvec(iicount)=1.0d0/ww
               enddo
            enddo
        enddo !n2 loop

!======================================================================
        ! Triple Integration
!======================================================================

        dlimn30 = 0.0D0
        dlimn31 = pi/4.0D0
        call gaussp(dlimn30,dlimn31,nn,xn3,wn3)
        do k=1,nn
           dlimn20 = atan(1.0D0/sin(xn3(k)))
           dlimn21 = pi/2.0D0
           call gaussp(dlimn20,dlimn21,nn,xn2,wn2)
           do j=1,nn
              dlimn10 = 0.0D0
              dlimn11 = 1.0D0/(alp(m)*sin(xn2(j))*cos(xn3(k)))
              call gaussp(dlimn10,dlimn11,nn,xn1,wn1)
              do i=1,nn
! Term 1
                 rr=xn1(i)
                 wttfact=(rr*rr*sin(xn2(j)))**(1.0d0/3.0d0)
  
                 uu=rr*sin(xn2(j))*cos(xn3(k))
                 vv=rr*sin(xn2(j))*sin(xn3(k))
                 ww=rr*cos(xn2(j))
  
                 iicount=iicount+1

                 wxvec(iicount)=(wn1(i)*wttfact)/(uu**2.0d0)
                 wyvec(iicount)=(wn2(j)*wttfact)/(vv**2.0d0)
                 wzvec(iicount)=(wn3(k)*wttfact)/(ww**2.0d0)
                 xvec(iicount)=1.0d0/uu
                 yvec(iicount)=1.0d0/vv
                 zvec(iicount)=1.0d0/ww

! Term 2
                 uu=rr*sin(xn2(j))*sin(xn3(k))
                 vv=rr*sin(xn2(j))*cos(xn3(k))
                 ww=rr*cos(xn2(j))
  
                 iicount=iicount+1

                 wxvec(iicount)=(wn1(i)*wttfact)/(uu**2.0d0)
                 wyvec(iicount)=(wn2(j)*wttfact)/(vv**2.0d0)
                 wzvec(iicount)=(wn3(k)*wttfact)/(ww**2.0d0)
                 xvec(iicount)=1.0d0/uu
                 yvec(iicount)=1.0d0/vv
                 zvec(iicount)=1.0d0/ww

! Term   3
                 uu=rr*sin(xn2(j))*sin(xn3(k))
                 vv=rr*cos(xn2(j))
                 ww=rr*sin(xn2(j))*cos(xn3(k))
   
                 iicount=iicount+1

                 wxvec(iicount)=(wn1(i)*wttfact)/(uu**2.0d0)
                 wyvec(iicount)=(wn2(j)*wttfact)/(vv**2.0d0)
                 wzvec(iicount)=(wn3(k)*wttfact)/(ww**2.0d0)
                 xvec(iicount)=1.0d0/uu
                 yvec(iicount)=1.0d0/vv
                 zvec(iicount)=1.0d0/ww

! Term 4
                 uu=rr*sin(xn2(j))*cos(xn3(k))
                 vv=rr*cos(xn2(j))
                 ww=rr*sin(xn2(j))*sin(xn3(k))

                 iicount=iicount+1

                 wxvec(iicount)=(wn1(i)*wttfact)/(uu**2.0d0)
                 wyvec(iicount)=(wn2(j)*wttfact)/(vv**2.0d0)
                 wzvec(iicount)=(wn3(k)*wttfact)/(ww**2.0d0)
                 xvec(iicount)=1.0d0/uu
                 yvec(iicount)=1.0d0/vv
                 zvec(iicount)=1.0d0/ww
  
! Term 5
                 xx=rr*cos(xn2(j))
                 vv=rr*sin(xn2(j))*cos(xn3(k))
                 ww=rr*sin(xn2(j))*sin(xn3(k))

                 iicount=iicount+1

                 wxvec(iicount)=(wn1(i)*wttfact)/(uu**2.0d0)
                 wyvec(iicount)=(wn2(j)*wttfact)/(vv**2.0d0)
                 wzvec(iicount)=(wn3(k)*wttfact)/(ww**2.0d0)
                 xvec(iicount)=1.0d0/uu
                 yvec(iicount)=1.0d0/vv
                 zvec(iicount)=1.0d0/ww
  
! Term 6
                 xx=rr*cos(xn2(j))
                 vv=rr*sin(xn2(j))*sin(xn3(k))
                 ww=rr*sin(xn2(j))*cos(xn3(k))
   
                 iicount=iicount+1

                 wxvec(iicount)=(wn1(i)*wttfact)/(uu**2.0d0)
                 wyvec(iicount)=(wn2(j)*wttfact)/(vv**2.0d0)
                 wzvec(iicount)=(wn3(k)*wttfact)/(ww**2.0d0)
                 xvec(iicount)=1.0d0/uu
                 yvec(iicount)=1.0d0/vv
                 zvec(iicount)=1.0d0/ww
          enddo
        enddo       
       enddo       

  505  continue ! m loop started at the top ends here 

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

      jjcount=0
      m=2
      do ii=1,iicount

        xnn=xvec(ii)
        ynn=yvec(ii)
        znn=zvec(ii)

!1------------- 1 0 0------------------------------
        if(xnn.ne.0.0d0.and.ynn.eq.0.0d0.and.znn.eq.0.0d0)then

          jjcount=jjcount+1

          xnn2=xnn
          ynn2=ynn
          znn2=znn

          wtt(jjcount,m)=wxvec(ii)*wyvec(ii)*wzvec(ii)
          rtvec(:)=(r0(:)-xnn2*tmat(1,:)-ynn2
     &         *tmat(2,:)-znn2*tmat(3,:))
          absrtvec=dsqrt(sum(rtvec(:)*rtvec(:)))

          do ix=1,3
             unitvec(ix,jjcount,m)=rtvec(ix)/absrtvec
          enddo
          denom(jjcount,m)=(sum(rtvec*rtvec))!**(dfloat(m)/2.0d0)
          !print 507,jjcount,xnn2,ynn2,znn2,absrtvec,wtt(jjcount,m)
!1------------- -1 0 0------------------------------
          jjcount=jjcount+1

          xnn2=-xnn
          ynn2=ynn
          znn2=znn

          wtt(jjcount,m)=wxvec(ii)*wyvec(ii)*wzvec(ii)
          rtvec(:)=(r0(:)-xnn2*tmat(1,:)-ynn2
     &         *tmat(2,:)-znn2*tmat(3,:))
          absrtvec=dsqrt(sum(rtvec(:)*rtvec(:)))

          do ix=1,3
             unitvec(ix,jjcount,m)=rtvec(ix)/absrtvec
          enddo
          denom(jjcount,m)=(sum(rtvec*rtvec))!**(dfloat(m)/2.0d0)
          !print 507,jjcount,xnn2,ynn2,znn2,absrtvec,wtt(jjcount,m)
        endif

!1------------- 0 1 0------------------------------
        if(xnn.eq.0.0d0.and.ynn.ne.0.0d0.and.znn.eq.0.0d0)then

          jjcount=jjcount+1

          xnn2=xnn
          ynn2=ynn
          znn2=znn

          wtt(jjcount,m)=wxvec(ii)*wyvec(ii)*wzvec(ii)
          rtvec(:)=(r0(:)-xnn2*tmat(1,:)-ynn2
     &         *tmat(2,:)-znn2*tmat(3,:))
          absrtvec=dsqrt(sum(rtvec(:)*rtvec(:)))

          do ix=1,3
             unitvec(ix,jjcount,m)=rtvec(ix)/absrtvec
          enddo
          denom(jjcount,m)=(sum(rtvec*rtvec))!**(dfloat(m)/2.0d0)
          !print 507,jjcount,xnn2,ynn2,znn2,absrtvec,wtt(jjcount,m)

!1------------- 0 -1 0------------------------------
          jjcount=jjcount+1

          xnn2=xnn
          ynn2=-ynn
          znn2=znn

          wtt(jjcount,m)=wxvec(ii)*wyvec(ii)*wzvec(ii)
          rtvec(:)=(r0(:)-xnn2*tmat(1,:)-ynn2
     &         *tmat(2,:)-znn2*tmat(3,:))
          absrtvec=dsqrt(sum(rtvec(:)*rtvec(:)))

          do ix=1,3
             unitvec(ix,jjcount,m)=rtvec(ix)/absrtvec
          enddo
          denom(jjcount,m)=(sum(rtvec*rtvec))!**(dfloat(m)/2.0d0)
          !print 507,jjcount,xnn2,ynn2,znn2,absrtvec,wtt(jjcount,m)
        endif

!1------------- 0 0 1------------------------------
        if(xnn.eq.0.0d0.and.ynn.eq.0.0d0.and.znn.ne.0.0d0)then

          jjcount=jjcount+1

          xnn2=xnn
          ynn2=ynn
          znn2=znn

          wtt(jjcount,m)=wxvec(ii)*wyvec(ii)*wzvec(ii)
          rtvec(:)=(r0(:)-xnn2*tmat(1,:)-ynn2
     &         *tmat(2,:)-znn2*tmat(3,:))
          absrtvec=dsqrt(sum(rtvec(:)*rtvec(:)))

          do ix=1,3
             unitvec(ix,jjcount,m)=rtvec(ix)/absrtvec
          enddo
          denom(jjcount,m)=(sum(rtvec*rtvec))!**(dfloat(m)/2.0d0)
          !print 507,jjcount,xnn2,ynn2,znn2,absrtvec,wtt(jjcount,m)

!1------------- 0 0 -1------------------------------
          jjcount=jjcount+1

          xnn2=xnn
          ynn2=ynn
          znn2=-znn

          wtt(jjcount,m)=wxvec(ii)*wyvec(ii)*wzvec(ii)
          rtvec(:)=(r0(:)-xnn2*tmat(1,:)-ynn2
     &         *tmat(2,:)-znn2*tmat(3,:))
          absrtvec=dsqrt(sum(rtvec(:)*rtvec(:)))

          do ix=1,3
             unitvec(ix,jjcount,m)=rtvec(ix)/absrtvec
          enddo
          denom(jjcount,m)=(sum(rtvec*rtvec))!**(dfloat(m)/2.0d0)
          !print 507,jjcount,xnn2,ynn2,znn2,absrtvec,wtt(jjcount,m)
        endif

!1------------- 0 1 1------------------------------
        if(xnn.eq.0.0d0.and.ynn.ne.0.0d0.and.znn.ne.0.0d0)then

          jjcount=jjcount+1

          xnn2=xnn
          ynn2=ynn
          znn2=znn

          wtt(jjcount,m)=wxvec(ii)*wyvec(ii)*wzvec(ii)
          rtvec(:)=(r0(:)-xnn2*tmat(1,:)-ynn2
     &         *tmat(2,:)-znn2*tmat(3,:))
          absrtvec=dsqrt(sum(rtvec(:)*rtvec(:)))

          do ix=1,3
             unitvec(ix,jjcount,m)=rtvec(ix)/absrtvec
          enddo
          denom(jjcount,m)=(sum(rtvec*rtvec))!**(dfloat(m)/2.0d0)
          !print 507,jjcount,xnn2,ynn2,znn2,absrtvec,wtt(jjcount,m)

!1------------- 0 -1 1------------------------------
          jjcount=jjcount+1

          xnn2=xnn
          ynn2=-ynn
          znn2=znn

          wtt(jjcount,m)=wxvec(ii)*wyvec(ii)*wzvec(ii)
          rtvec(:)=(r0(:)-xnn2*tmat(1,:)-ynn2
     &         *tmat(2,:)-znn2*tmat(3,:))
          absrtvec=dsqrt(sum(rtvec(:)*rtvec(:)))

          do ix=1,3
             unitvec(ix,jjcount,m)=rtvec(ix)/absrtvec
          enddo
          denom(jjcount,m)=(sum(rtvec*rtvec))!**(dfloat(m)/2.0d0)
          !print 507,jjcount,xnn2,ynn2,znn2,absrtvec,wtt(jjcount,m)

!1------------- 0 1 -1------------------------------
          jjcount=jjcount+1

          xnn2=xnn
          ynn2=ynn
          znn2=-znn

          wtt(jjcount,m)=wxvec(ii)*wyvec(ii)*wzvec(ii)
          rtvec(:)=(r0(:)-xnn2*tmat(1,:)-ynn2
     &         *tmat(2,:)-znn2*tmat(3,:))
          absrtvec=dsqrt(sum(rtvec(:)*rtvec(:)))

          do ix=1,3
             unitvec(ix,jjcount,m)=rtvec(ix)/absrtvec
          enddo
          denom(jjcount,m)=(sum(rtvec*rtvec))!**(dfloat(m)/2.0d0)
          !print 507,jjcount,xnn2,ynn2,znn2,absrtvec,wtt(jjcount,m)

!1------------- 0 -1 -1------------------------------
          jjcount=jjcount+1

          xnn2=xnn
          ynn2=-ynn
          znn2=-znn

          wtt(jjcount,m)=wxvec(ii)*wyvec(ii)*wzvec(ii)
          rtvec(:)=(r0(:)-xnn2*tmat(1,:)-ynn2
     &         *tmat(2,:)-znn2*tmat(3,:))
          absrtvec=dsqrt(sum(rtvec(:)*rtvec(:)))

          do ix=1,3
             unitvec(ix,jjcount,m)=rtvec(ix)/absrtvec
          enddo
          denom(jjcount,m)=(sum(rtvec*rtvec))!**(dfloat(m)/2.0d0)
          !print 507,jjcount,xnn2,ynn2,znn2,absrtvec,wtt(jjcount,m)

        endif
!1------------- 1 0 1------------------------------
        if(xnn.ne.0.0d0.and.ynn.eq.0.0d0.and.znn.ne.0.0d0)then

          jjcount=jjcount+1

          xnn2=xnn
          ynn2=ynn
          znn2=znn

          wtt(jjcount,m)=wxvec(ii)*wyvec(ii)*wzvec(ii)
          rtvec(:)=(r0(:)-xnn2*tmat(1,:)-ynn2
     &         *tmat(2,:)-znn2*tmat(3,:))
          absrtvec=dsqrt(sum(rtvec(:)*rtvec(:)))

          do ix=1,3
             unitvec(ix,jjcount,m)=rtvec(ix)/absrtvec
          enddo
          denom(jjcount,m)=(sum(rtvec*rtvec))!**(dfloat(m)/2.0d0)
          !print 507,jjcount,xnn2,ynn2,znn2,absrtvec,wtt(jjcount,m)

!1------------- -1 0 1------------------------------
          jjcount=jjcount+1

          xnn2=-xnn
          ynn2=ynn
          znn2=znn

          wtt(jjcount,m)=wxvec(ii)*wyvec(ii)*wzvec(ii)
          rtvec(:)=(r0(:)-xnn2*tmat(1,:)-ynn2
     &         *tmat(2,:)-znn2*tmat(3,:))
          absrtvec=dsqrt(sum(rtvec(:)*rtvec(:)))

          do ix=1,3
             unitvec(ix,jjcount,m)=rtvec(ix)/absrtvec
          enddo
          denom(jjcount,m)=(sum(rtvec*rtvec))!**(dfloat(m)/2.0d0)
          !print 507,jjcount,xnn2,ynn2,znn2,absrtvec,wtt(jjcount,m)

!1------------- 1 0 -1------------------------------
          jjcount=jjcount+1

          xnn2=xnn
          ynn2=ynn
          znn2=-znn

          wtt(jjcount,m)=wxvec(ii)*wyvec(ii)*wzvec(ii)
          rtvec(:)=(r0(:)-xnn2*tmat(1,:)-ynn2
     &         *tmat(2,:)-znn2*tmat(3,:))
          absrtvec=dsqrt(sum(rtvec(:)*rtvec(:)))

          do ix=1,3
             unitvec(ix,jjcount,m)=rtvec(ix)/absrtvec
          enddo
          denom(jjcount,m)=(sum(rtvec*rtvec))!**(dfloat(m)/2.0d0)
          !print 507,jjcount,xnn2,ynn2,znn2,absrtvec,wtt(jjcount,m)

!1------------- -1 0 -1------------------------------
          jjcount=jjcount+1

          xnn2=-xnn
          ynn2=ynn
          znn2=-znn

          wtt(jjcount,m)=wxvec(ii)*wyvec(ii)*wzvec(ii)
          rtvec(:)=(r0(:)-xnn2*tmat(1,:)-ynn2
     &         *tmat(2,:)-znn2*tmat(3,:))
          absrtvec=dsqrt(sum(rtvec(:)*rtvec(:)))

          do ix=1,3
             unitvec(ix,jjcount,m)=rtvec(ix)/absrtvec
          enddo
          denom(jjcount,m)=(sum(rtvec*rtvec))!**(dfloat(m)/2.0d0)
          !print 507,jjcount,xnn2,ynn2,znn2,absrtvec,wtt(jjcount,m)

        endif
!1------------- 1 1 0------------------------------

        if(xnn.ne.0.0d0.and.ynn.ne.0.0d0.and.znn.eq.0.0d0)then

          jjcount=jjcount+1

          xnn2=xnn
          ynn2=ynn
          znn2=znn

          wtt(jjcount,m)=wxvec(ii)*wyvec(ii)*wzvec(ii)
          rtvec(:)=(r0(:)-xnn2*tmat(1,:)-ynn2
     &         *tmat(2,:)-znn2*tmat(3,:))
          absrtvec=dsqrt(sum(rtvec(:)*rtvec(:)))

          do ix=1,3
             unitvec(ix,jjcount,m)=rtvec(ix)/absrtvec
          enddo
          denom(jjcount,m)=(sum(rtvec*rtvec))!**(dfloat(m)/2.0d0)
          !print 507,jjcount,xnn2,ynn2,znn2,absrtvec,wtt(jjcount,m)

!1------------- -1 1 0------------------------------
          jjcount=jjcount+1

          xnn2=-xnn
          ynn2=ynn
          znn2=znn

          wtt(jjcount,m)=wxvec(ii)*wyvec(ii)*wzvec(ii)
          rtvec(:)=(r0(:)-xnn2*tmat(1,:)-ynn2
     &         *tmat(2,:)-znn2*tmat(3,:))
          absrtvec=dsqrt(sum(rtvec(:)*rtvec(:)))

          do ix=1,3
             unitvec(ix,jjcount,m)=rtvec(ix)/absrtvec
          enddo
          denom(jjcount,m)=(sum(rtvec*rtvec))!**(dfloat(m)/2.0d0)
          !print 507,jjcount,xnn2,ynn2,znn2,absrtvec,wtt(jjcount,m)

!1------------- 1 -1 0------------------------------
          jjcount=jjcount+1

          xnn2=xnn
          ynn2=-ynn
          znn2=znn

          wtt(jjcount,m)=wxvec(ii)*wyvec(ii)*wzvec(ii)
          rtvec(:)=(r0(:)-xnn2*tmat(1,:)-ynn2
     &         *tmat(2,:)-znn2*tmat(3,:))
          absrtvec=dsqrt(sum(rtvec(:)*rtvec(:)))

          do ix=1,3
             unitvec(ix,jjcount,m)=rtvec(ix)/absrtvec
          enddo
          denom(jjcount,m)=(sum(rtvec*rtvec))!**(dfloat(m)/2.0d0)
          !print 507,jjcount,xnn2,ynn2,znn2,absrtvec,wtt(jjcount,m)

!1------------- -1 -1 0------------------------------
          jjcount=jjcount+1

          xnn2=-xnn
          ynn2=-ynn
          znn2=znn

          wtt(jjcount,m)=wxvec(ii)*wyvec(ii)*wzvec(ii)
          rtvec(:)=(r0(:)-xnn2*tmat(1,:)-ynn2
     &         *tmat(2,:)-znn2*tmat(3,:))
          absrtvec=dsqrt(sum(rtvec(:)*rtvec(:)))

          do ix=1,3
             unitvec(ix,jjcount,m)=rtvec(ix)/absrtvec
          enddo
          denom(jjcount,m)=(sum(rtvec*rtvec))!**(dfloat(m)/2.0d0)
          !print 507,jjcount,xnn2,ynn2,znn2,absrtvec,wtt(jjcount,m)
        endif

        if(xnn.ne.0.0d0.and.ynn.ne.0.0d0.and.znn.ne.0.0d0)then
!1------------- 1 1 1------------------------------
          jjcount=jjcount+1

          xnn2=xnn
          ynn2=ynn
          znn2=znn

          wtt(jjcount,m)=wxvec(ii)*wyvec(ii)*wzvec(ii)

          rtvec(:)=(r0(:)-xnn2*tmat(1,:)-ynn2
     &           *tmat(2,:)-znn2*tmat(3,:))
          absrtvec=dsqrt(sum(rtvec(:)*rtvec(:)))
          do ix=1,3
             unitvec(ix,jjcount,m)=rtvec(ix)/absrtvec
          enddo
          denom(jjcount,m)=(sum(rtvec*rtvec))!**(dfloat(m)/2.0d0)

          !print 507,jjcount,xnn2,ynn2,znn2,absrtvec,wtt(jjcount,m)
!2------------- 1 1 -1------------------------------

          jjcount=jjcount+1

          xnn2=xnn
          ynn2=ynn
          znn2=-znn

          wtt(jjcount,m)=wxvec(ii)*wyvec(ii)*wzvec(ii)

          rtvec(:)=(r0(:)-xnn2*tmat(1,:)-ynn2
     &           *tmat(2,:)-znn2*tmat(3,:))
          absrtvec=dsqrt(sum(rtvec(:)*rtvec(:)))
          do ix=1,3
             unitvec(ix,jjcount,m)=rtvec(ix)/absrtvec
          enddo
          denom(jjcount,m)=(sum(rtvec*rtvec))!**(dfloat(m)/2.0d0)
          !print 507,jjcount,xnn2,ynn2,znn2,absrtvec,wtt(jjcount,m)

!3------------- 1 -1 1------------------------------

          jjcount=jjcount+1

          xnn2=xnn
          ynn2=-ynn
          znn2=znn

          wtt(jjcount,m)=wxvec(ii)*wyvec(ii)*wzvec(ii)

          rtvec(:)=(r0(:)-xnn2*tmat(1,:)-ynn2
     &           *tmat(2,:)-znn2*tmat(3,:))
          absrtvec=dsqrt(sum(rtvec(:)*rtvec(:)))
          do ix=1,3
             unitvec(ix,jjcount,m)=rtvec(ix)/absrtvec
          enddo
          denom(jjcount,m)=(sum(rtvec*rtvec))!**(dfloat(m)/2.0d0)
          !print 507,jjcount,xnn2,ynn2,znn2,absrtvec,wtt(jjcount,m)

!4------------- -1 1 1------------------------------

          jjcount=jjcount+1

          xnn2=-xnn
          ynn2=ynn
          znn2=znn

          wtt(jjcount,m)=wxvec(ii)*wyvec(ii)*wzvec(ii)

          rtvec(:)=(r0(:)-xnn2*tmat(1,:)-ynn2
     &           *tmat(2,:)-znn2*tmat(3,:))
          absrtvec=dsqrt(sum(rtvec(:)*rtvec(:)))
          do ix=1,3
             unitvec(ix,jjcount,m)=rtvec(ix)/absrtvec
          enddo
          denom(jjcount,m)=(sum(rtvec*rtvec))!**(dfloat(m)/2.0d0)
          !print 507,jjcount,xnn2,ynn2,znn2,absrtvec,wtt(jjcount,m)

!5------------- -1 -1 1------------------------------

          jjcount=jjcount+1

          xnn2=-xnn
          ynn2=-ynn
          znn2=znn

          wtt(jjcount,m)=wxvec(ii)*wyvec(ii)*wzvec(ii)

          rtvec(:)=(r0(:)-xnn2*tmat(1,:)-ynn2
     &           *tmat(2,:)-znn2*tmat(3,:))
          absrtvec=dsqrt(sum(rtvec(:)*rtvec(:)))
          do ix=1,3
             unitvec(ix,jjcount,m)=rtvec(ix)/absrtvec
          enddo
          denom(jjcount,m)=(sum(rtvec*rtvec))!**(dfloat(m)/2.0d0)
          !print 507,jjcount,xnn2,ynn2,znn2,absrtvec,wtt(jjcount,m)

!6------------- -1 1 -1------------------------------

          jjcount=jjcount+1

          xnn2=-xnn
          ynn2=ynn
          znn2=-znn

          wtt(jjcount,m)=wxvec(ii)*wyvec(ii)*wzvec(ii)

          rtvec(:)=(r0(:)-xnn2*tmat(1,:)-ynn2
     &           *tmat(2,:)-znn2*tmat(3,:))
          absrtvec=dsqrt(sum(rtvec(:)*rtvec(:)))
          do ix=1,3
             unitvec(ix,jjcount,m)=rtvec(ix)/absrtvec
          enddo
          denom(jjcount,m)=(sum(rtvec*rtvec))!**(dfloat(m)/2.0d0)
          !print 507,jjcount,xnn2,ynn2,znn2,absrtvec,wtt(jjcount,m)

!7-------------1 -1 -1------------------------------

          jjcount=jjcount+1

          xnn2=xnn
          ynn2=-ynn
          znn2=-znn

          wtt(jjcount,m)=wxvec(ii)*wyvec(ii)*wzvec(ii)

          rtvec(:)=(r0(:)-xnn2*tmat(1,:)-ynn2
     &           *tmat(2,:)-znn2*tmat(3,:))
          absrtvec=dsqrt(sum(rtvec(:)*rtvec(:)))
          do ix=1,3
             unitvec(ix,jjcount,m)=rtvec(ix)/absrtvec
          enddo
          denom(jjcount,m)=(sum(rtvec*rtvec))!**(dfloat(m)/2.0d0)
          !print 507,jjcount,xnn2,ynn2,znn2,absrtvec,wtt(jjcount,m)

!8------------- -1 -1 -1------------------------------

          jjcount=jjcount+1

          xnn2=-xnn
          ynn2=-ynn
          znn2=-znn

          wtt(jjcount,m)=wxvec(ii)*wyvec(ii)*wzvec(ii)

          rtvec(:)=(r0(:)-xnn2*tmat(1,:)-ynn2
     &           *tmat(2,:)-znn2*tmat(3,:))
          absrtvec=dsqrt(sum(rtvec(:)*rtvec(:)))
          do ix=1,3
             unitvec(ix,jjcount,m)=rtvec(ix)/absrtvec
          enddo
          denom(jjcount,m)=(sum(rtvec*rtvec))!**(dfloat(m)/2.0d0)
          !print 507,jjcount,xnn2,ynn2,znn2,absrtvec,wtt(jjcount,m)

        endif
        enddo

!      do jj=1,jjcount
!        print *,(unitvec(ix,jj,m),ix=1,3)
!      enddo

!      if(igrid.eq.iicount)then
!        do m=2,itotmultipole
        ALLOCATE(unitvec2(3,jjcount))

        do jj=1,jjcount
        unitvec2(1,jj)=unitvec(1,jj,m)
        unitvec2(2,jj)=unitvec(2,jj,m)
        unitvec2(3,jj)=unitvec(3,jj,m)
        enddo

!        do jj=1,jjcount
!          print *,(unitvec2(ix,jj),ix=1,3)
!        enddo

        ALLOCATE(POLY(jjcount,(LM+1)**2))
        CALL HARMONICS(jjcount,jjcount,LM,unitvec2,POLY,NPOLY)

! Checking the output YLM generated above
!      print *,NPOLY
!      do i=1,NPOLY
!        do j=1,iicount
!          print *,POLY(j,i)
!        enddo
!      enddo
!
!      stop

      sumdirect=0.0d0

      do jj=1,jjcount
         do kk=1,ILMTOT
           lmm=ILM(kk)
           mm=ceiling(dsqrt(dfloat(lmm)))
           sumdirect=sumdirect+(SMULTI(kk)*POLY(jj,lmm)*wtt(jj,m))
     &             /(denom(jj,m))**(dfloat(mm)/2.0d0)
         enddo
      enddo
      !print *,'sum',jjcount,sumdirect
      print *,'Sum+Int=',sumdirect, jjcount

      DEALLOCATE(POLY,unitvec,denom,ILM,wtt)

      close(399)
      close(299)

 507   format(i3,6f14.4)
 10    format(i10,1x,15G15.10)
 11    format(i5,i5,1x,15G15.10)
      END SUBROUTINE LAT_SUM

!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!C
!C GAUSSP CREATES A N-POINT GAUSS-LEGENDRE MESH SUITABLE FOR THE 
!C INTEGRATION OF POLYNOMIALS OF DEGREE 2*N-1 IN THE INTERVAL (Y1,Y2)
!C
        SUBROUTINE GAUSSP(Y1,Y2,N,Y,WY)
        IMPLICIT REAL*8 (A-H,O-Z)
        DIMENSION Y(N),WY(N)
        PARAMETER (EPS=1.0D-14)
        PARAMETER (MAX=1000)
        DIMENSION X(MAX),W(MAX)
        SAVE
        DATA N_SAV/-1/
        DATA PI/3.14159265358979323844D0/
!C
        IF (N. LE. 0) RETURN
        IF (N .GT. MAX) THEN
         PRINT*,'GAUSSP: MAX MUST BE AT LEAST: ',N
         STOP          
         ! CALL STOPIT Changed to STOP by YS here
        END IF
        IF (N .EQ. 1) THEN
         WY(1)=Y2-Y1
         Y(1)=0.5D0*(Y2+Y1)
         RETURN
        END IF
!C
!C FIRST TRY PRESTORED MESHES
!C IF NONE OF THEM MATCHES, DEFINE MESH VIA ITERATION OF CHEBYSHEV POINTS
!C POINTS/WEIGHTS WILL BE DEFINED FOR INTERVAL (-1:1)
!C
        IF (N .NE. N_SAV) THEN
         N_SAV=N
         M=N
         CALL GAUSSPP(M,X,W)
         IF (M .NE. N) THEN
          M=(N+1)/2
          DO I=1,M
           Y(I)=COS(PI*(I-0.25D0)/(N+0.5D0))
          END DO
          DO 100 I=1,M
           Z=Y(I)
  20       CONTINUE
            P1=1.0D0
            P2=0.0D0
            DO J=1,N
             P3=P2
             P2=P1
             P1=((2*J-1)*Z*P2-(J-1)*P3)/J
            END DO
            PP=N*(Z*P1-P2)/(Z*Z-1.0D0)
            Z1=Z
            Z=Z1-P1/PP
            IF (ABS(Z-Z1).GT.EPS) GOTO 20
           CONTINUE
           X(I)    = -Z
           X(N+1-I)=  Z
           W(I)=2.0D0/((1.0D0-Z*Z)*PP*PP)
           W(N+1-I)=W(I)
 100      CONTINUE
         END IF
        END IF
!C
!C SCALE WEIGHTS AND POINTS:
!C
        DO I=1,N
         FACT= 0.5D0*(Y2-Y1)
         Y (I)= Y1+FACT*(X(I)+1.0D0)
         WY(I)= FACT*W(I)
        END DO
        RETURN
        END
!C
!C ******************************************************************
!C


      SUBROUTINE GAUSSPP(NPT,XPT,WHT)
      IMPLICIT REAL*8 (A-H,O-Z)
      DIMENSION XPT(*),WHT(*)
!C
!C     THIS ROUTINE SETS UP ABSCISSAE AND WEIGHTS FOR NPT-POINT
!C       GAUSS-LEGENDRE INTEGRATION IN THE INTERVAL (-1,1).
!C
!C     ON RETURN, THE FUNCTION TO BE INTEGRATED SHOULD BE EVALUATED
!C       AT THE POINTS XPT(I).  INTEGRAL = SUM(I=1,NPT) F(XPT(I))*WHT(I).
!C
!C     FOR NPT=1-14:  VALUES FROM Z. KOPAL, NUMERICAL ANALYSIS, 1961, A4
!C         NPT=16-64: FROM ALCHEMY GAUSS ROUTINE
!C         NPT=96:    FROM DECK OF R. NERF JUNE 1973.
!C
      DIMENSION IOK(22),IST(22)
      DIMENSION X(239),W(239)
      DIMENSION X2(1),W2(1)
      DIMENSION X3(2),W3(2)
      DIMENSION X4(2),W4(2)
      DIMENSION X5(3),W5(3)
      DIMENSION X6(3),W6(3)
      DIMENSION X7(4),W7(4)
      DIMENSION X8(4),W8(4)
      DIMENSION X9(5),W9(5)
      DIMENSION X10(5),W10(5)
      DIMENSION X11(6),W11(6)
      DIMENSION X12(6),W12(6)
      DIMENSION X13(7),W13(7)
      DIMENSION X14(7),W14(7)
      DIMENSION X16(8),W16(8)
      DIMENSION X20(10),W20(10)
      DIMENSION X24(12),W24(12)
      DIMENSION X28(14),W28(14)
      DIMENSION X32(16),W32(16)
      DIMENSION X40(20),W40(20)
      DIMENSION X48(24),W48(24)
      DIMENSION X64(32),W64(32)
      DIMENSION X96(48),W96(48)
      EQUIVALENCE(X2(1), X(  1)),(W2(1), W(  1))
      EQUIVALENCE(X3(1), X(  2)),(W3(1), W(  2))
      EQUIVALENCE(X4(1), X(  4)),(W4(1), W(  4))
      EQUIVALENCE(X5(1), X(  6)),(W5(1), W(  6))
      EQUIVALENCE(X6(1), X(  9)),(W6(1), W(  9))
      EQUIVALENCE(X7(1), X( 12)),(W7(1), W( 12))
      EQUIVALENCE(X8(1), X( 16)),(W8(1), W( 16))
      EQUIVALENCE(X9(1), X( 20)),(W9(1), W( 20))
      EQUIVALENCE(X10(1),X( 25)),(W10(1),W( 25))
      EQUIVALENCE(X11(1),X( 30)),(W11(1),W( 30))
      EQUIVALENCE(X12(1),X( 36)),(W12(1),W( 36))
      EQUIVALENCE(X13(1),X( 42)),(W13(1),W( 42))
      EQUIVALENCE(X14(1),X( 49)),(W14(1),W( 49))
      EQUIVALENCE(X16(1),X( 56)),(W16(1),W( 56))
      EQUIVALENCE(X20(1),X( 64)),(W20(1),W( 64))
      EQUIVALENCE(X24(1),X( 74)),(W24(1),W( 74))
      EQUIVALENCE(X28(1),X( 86)),(W28(1),W( 86))
      EQUIVALENCE(X32(1),X(100)),(W32(1),W(100))
      EQUIVALENCE(X40(1),X(116)),(W40(1),W(116))
      EQUIVALENCE(X48(1),X(136)),(W48(1),W(136))
      EQUIVALENCE(X64(1),X(160)),(W64(1),W(160))
      EQUIVALENCE(X96(1),X(192)),(W96(1),W(192))
      SAVE
!C
!C MX: TOTAL NUMBER OF STORED MESHES
!C IOK: NUMBER OF POINTS FOR STORED MESHES
!C IST: START INDEX FOR STORED MESHES IN ARRAYS X,W
!C
      DATA MX/22/
      DATA IOK/2,3,4,5,6,7,8,9,10,11,12,13,
     &         14,16,20,24,28,32,40,48,64,96/
      DATA IST/1,2,4,6,9,12,16,20,25,30,36,42,
     &         49,56,64,74,86,100,116,136,160,192/
!C
      DATA X2 ( 1)/0.577350269189626D+0/,W2( 1) /0.999999999999999D+0/
      DATA X3 ( 1)/0.774596669241483D+0/,W3( 1) /0.555555555555556D+0/
      DATA X3 ( 2)/0.0              D+0/,W3( 2) /0.888888888888889D+0/
      DATA X4 ( 1)/0.861136311594053D+0/,W4( 1) /0.347854845137454D+0/
      DATA X4 ( 2)/0.339981043584856D+0/,W4( 2) /0.652145154862546D+0/
      DATA X5 ( 1)/0.906179845938664D+0/,W5( 1) /0.236926885056189D+0/
      DATA X5 ( 2)/0.538469310105683D+0/,W5( 2) /0.478628670499366D+0/
      DATA X5 ( 3)/0.0              D+0/,W5( 3) /0.568888888888889D+0/
      DATA X6 ( 1)/0.932469514203152D+0/,W6( 1) /0.171324492379170D+0/
      DATA X6 ( 2)/0.661209386466265D+0/,W6( 2) /0.360761573048139D+0/
      DATA X6 ( 3)/0.238619186083197D+0/,W6( 3) /0.467913934572691D+0/
      DATA X7 ( 1)/0.949107912342759D+0/,W7( 1) /0.129484966168870D+0/
      DATA X7 ( 2)/0.741531185599394D+0/,W7( 2) /0.279705391489277D+0/
      DATA X7 ( 3)/0.405845151377397D+0/,W7( 3) /0.381830050505119D+0/
      DATA X7 ( 4)/0.0              D+0/,W7( 4) /0.417959183673469D+0/
      DATA X8 ( 1)/0.960289856497536D+0/,W8( 1) /0.101228536290376D+0/
      DATA X8 ( 2)/0.796666477413627D+0/,W8( 2) /0.222381034453374D+0/
      DATA X8 ( 3)/0.525532409916329D+0/,W8( 3) /0.313706645877887D+0/
      DATA X8 ( 4)/0.183434642495650D+0/,W8( 4) /0.362683783378362D+0/
      DATA X9 ( 1)/0.968160239507626D+0/,W9( 1) /0.812743883615739D-1/
      DATA X9 ( 2)/0.836031107326636D+0/,W9( 2) /0.180648160694857D+0/
      DATA X9 ( 3)/0.613371432700590D+0/,W9( 3) /0.260610696402935D+0/
      DATA X9 ( 4)/0.324253423403809D+0/,W9( 4) /0.312347077040003D+0/
      DATA X9 ( 5)/0.0              D+0/,W9( 5) /0.330239355001260D+0/
      DATA X10( 1)/0.973906528517172D+0/,W10( 1)/0.666713443086879D-1/
      DATA X10( 2)/0.865063366688985D+0/,W10( 2)/0.149451349150581D+0/
      DATA X10( 3)/0.679409568299024D+0/,W10( 3)/0.219086362515982D+0/
      DATA X10( 4)/0.433395394129247D+0/,W10( 4)/0.269266719309996D+0/
      DATA X10( 5)/0.148874338981631D+0/,W10( 5)/0.295524224714753D+0/
      DATA X11( 1)/0.978228658146057D+0/,W11( 1)/0.556685671161740D-1/
      DATA X11( 2)/0.887062599768095D+0/,W11( 2)/0.125580369464905D+0/
      DATA X11( 3)/0.730152005574049D+0/,W11( 3)/0.186290210927734D+0/
      DATA X11( 4)/0.519096129206812D+0/,W11( 4)/0.233193764591990D+0/
      DATA X11( 5)/0.269543155952345D+0/,W11( 5)/0.262804544510247D+0/
      DATA X11( 6)/0.0              D+0/,W11( 6)/0.272925086777901D+0/
      DATA X12( 1)/0.981560634246719D+0/,W12( 1)/0.471753363865120D-1/
      DATA X12( 2)/0.904117256370475D+0/,W12( 2)/0.106939325995318D+0/
      DATA X12( 3)/0.769902674194305D+0/,W12( 3)/0.160078328543346D+0/
      DATA X12( 4)/0.587317954286617D+0/,W12( 4)/0.203167426723066D+0/
      DATA X12( 5)/0.367831498998180D+0/,W12( 5)/0.233492536538355D+0/
      DATA X12( 6)/0.125233408511469D+0/,W12( 6)/0.249147045813403D+0/
      DATA X13( 1)/0.984183054718588D+0/,W13( 1)/0.404840047653160D-1/
      DATA X13( 2)/0.917598399222978D+0/,W13( 2)/0.921214998377279D-1/
      DATA X13( 3)/0.801578090733310D+0/,W13( 3)/0.138873510219787D+0/
      DATA X13( 4)/0.642349339440340D+0/,W13( 4)/0.178145980761946D+0/
      DATA X13( 5)/0.448492751036447D+0/,W13( 5)/0.207816047536889D+0/
      DATA X13( 6)/0.230458315955135D+0/,W13( 6)/0.226283180262897D+0/
      DATA X13( 7)/0.0              D+0/,W13( 7)/0.232551553230874D+0/
      DATA X14( 1)/0.986283808696812D+0/,W14( 1)/0.351194603317520D-1/
      DATA X14( 2)/0.928434883663574D+0/,W14( 2)/0.801580871597599D-1/
      DATA X14( 3)/0.827201315069765D+0/,W14( 3)/0.121518570687903D+0/
      DATA X14( 4)/0.687292904811685D+0/,W14( 4)/0.157203167158194D+0/
      DATA X14( 5)/0.515248636358154D+0/,W14( 5)/0.185538397477938D+0/
      DATA X14( 6)/0.319112368927890D+0/,W14( 6)/0.205198463721296D+0/
      DATA X14( 7)/0.108054948707344D+0/,W14( 7)/0.215263853463158D+0/
      DATA X16( 1)/0.989400934991650D+0/,W16( 1)/0.271524594117540D-1/
      DATA X16( 2)/0.944575023073232D+0/,W16( 2)/0.622535239386480D-1/
      DATA X16( 3)/0.865631202387832D+0/,W16( 3)/0.951585116824929D-1/
      DATA X16( 4)/0.755404408355003D+0/,W16( 4)/0.124628971255534D+0/
      DATA X16( 5)/0.617876244402644D+0/,W16( 5)/0.149595988816577D+0/
      DATA X16( 6)/0.458016777657227D+0/,W16( 6)/0.169156519395002D+0/
      DATA X16( 7)/0.281603550779259D+0/,W16( 7)/0.182603415044923D+0/
      DATA X16( 8)/0.950125098376369D-1/,W16( 8)/0.189450610455068D+0/
      DATA X20( 1)/0.993128599185095D+0/,W20( 1)/0.176140071391520D-1/
      DATA X20( 2)/0.963971927277914D+0/,W20( 2)/0.406014298003870D-1/
      DATA X20( 3)/0.912234428251326D+0/,W20( 3)/0.626720483341089D-1/
      DATA X20( 4)/0.839116971822219D+0/,W20( 4)/0.832767415767049D-1/
      DATA X20( 5)/0.746331906460151D+0/,W20( 5)/0.101930119817240D+0/
      DATA X20( 6)/0.636053680726515D+0/,W20( 6)/0.118194531961518D+0/
      DATA X20( 7)/0.510867001950827D+0/,W20( 7)/0.131688638449177D+0/
      DATA X20( 8)/0.373706088715419D+0/,W20( 8)/0.142096109318382D+0/
      DATA X20( 9)/0.227785851141645D+0/,W20( 9)/0.149172986472604D+0/
      DATA X20(10)/0.765265211334969D-1/,W20(10)/0.152753387130726D+0/
      DATA X24( 1)/0.995187219997021D+0/,W24( 1)/0.123412297999870D-1/
      DATA X24( 2)/0.974728555971309D+0/,W24( 2)/0.285313886289340D-1/
      DATA X24( 3)/0.938274552002733D+0/,W24( 3)/0.442774388174200D-1/
      DATA X24( 4)/0.886415527004401D+0/,W24( 4)/0.592985849154370D-1/
      DATA X24( 5)/0.820001985973903D+0/,W24( 5)/0.733464814110799D-1/
      DATA X24( 6)/0.740124191578554D+0/,W24( 6)/0.861901615319529D-1/
      DATA X24( 7)/0.648093651936975D+0/,W24( 7)/0.976186521041139D-1/
      DATA X24( 8)/0.545421471388839D+0/,W24( 8)/0.107444270115966D+0/
      DATA X24( 9)/0.433793507626045D+0/,W24( 9)/0.115505668053726D+0/
      DATA X24(10)/0.315042679696163D+0/,W24(10)/0.121670472927803D+0/
      DATA X24(11)/0.191118867473616D+0/,W24(11)/0.125837456346828D+0/
      DATA X24(12)/0.640568928626059D-1/,W24(12)/0.127938195346752D+0/
      DATA X28( 1)/0.996442497573954D+0/,W28( 1)/0.912428259309400D-2/
      DATA X28( 2)/0.981303165370873D+0/,W28( 2)/0.211321125927710D-1/
      DATA X28( 3)/0.954259280628938D+0/,W28( 3)/0.329014277823040D-1/
      DATA X28( 4)/0.915633026392132D+0/,W28( 4)/0.442729347590040D-1/
      DATA X28( 5)/0.865892522574395D+0/,W28( 5)/0.551073456757170D-1/
      DATA X28( 6)/0.805641370917179D+0/,W28( 6)/0.652729239669989D-1/
      DATA X28( 7)/0.735610878013632D+0/,W28( 7)/0.746462142345689D-1/
      DATA X28( 8)/0.656651094038865D+0/,W28( 8)/0.831134172289009D-1/
      DATA X28( 9)/0.569720471811402D+0/,W28( 9)/0.905717443930329D-1/
      DATA X28(10)/0.475874224955118D+0/,W28(10)/0.969306579979299D-1/
      DATA X28(11)/0.376251516089079D+0/,W28(11)/0.102112967578061D+0/
      DATA X28(12)/0.272061627635178D+0/,W28(12)/0.106055765922846D+0/
      DATA X28(13)/0.164569282133381D+0/,W28(13)/0.108711192258294D+0/
      DATA X28(14)/0.550792898840340D-1/,W28(14)/0.110047013016475D+0/
      DATA X32( 1)/0.997263861849481D+0/,W32( 1)/0.701861000947000D-2/
      DATA X32( 2)/0.985611511545268D+0/,W32( 2)/0.162743947309060D-1/
      DATA X32( 3)/0.964762255587506D+0/,W32( 3)/0.253920653092620D-1/
      DATA X32( 4)/0.934906075937740D+0/,W32( 4)/0.342738629130210D-1/
      DATA X32( 5)/0.896321155766052D+0/,W32( 5)/0.428358980222270D-1/
      DATA X32( 6)/0.849367613732570D+0/,W32( 6)/0.509980592623760D-1/
      DATA X32( 7)/0.794483795967942D+0/,W32( 7)/0.586840934785350D-1/
      DATA X32( 8)/0.732182118740290D+0/,W32( 8)/0.658222227763619D-1/
      DATA X32( 9)/0.663044266930215D+0/,W32( 9)/0.723457941088479D-1/
      DATA X32(10)/0.587715757240762D+0/,W32(10)/0.781938957870699D-1/
      DATA X32(11)/0.506899908932229D+0/,W32(11)/0.833119242269469D-1/
      DATA X32(12)/0.421351276130635D+0/,W32(12)/0.876520930044039D-1/
      DATA X32(13)/0.331868602282128D+0/,W32(13)/0.911738786957639D-1/
      DATA X32(14)/0.239287362252137D+0/,W32(14)/0.938443990808039D-1/
      DATA X32(15)/0.144471961582796D+0/,W32(15)/0.956387200792749D-1/
      DATA X32(16)/0.483076656877380D-1/,W32(16)/0.965400885147279D-1/
      DATA X40( 1)/0.998237709710559D+0/,W40( 1)/0.452127709853300D-2/
      DATA X40( 2)/0.990726238699457D+0/,W40( 2)/0.104982845311530D-1/
      DATA X40( 3)/0.977259949983774D+0/,W40( 3)/0.164210583819080D-1/
      DATA X40( 4)/0.957916819213792D+0/,W40( 4)/0.222458491941670D-1/
      DATA X40( 5)/0.932812808278676D+0/,W40( 5)/0.279370069800230D-1/
      DATA X40( 6)/0.902098806968874D+0/,W40( 6)/0.334601952825480D-1/
      DATA X40( 7)/0.865959503212259D+0/,W40( 7)/0.387821679744720D-1/
      DATA X40( 8)/0.824612230833312D+0/,W40( 8)/0.438709081856730D-1/
      DATA X40( 9)/0.778305651426519D+0/,W40( 9)/0.486958076350720D-1/
      DATA X40(10)/0.727318255189927D+0/,W40(10)/0.532278469839370D-1/
      DATA X40(11)/0.671956684614179D+0/,W40(11)/0.574397690993910D-1/
      DATA X40(12)/0.612553889667980D+0/,W40(12)/0.613062424929290D-1/
      DATA X40(13)/0.549467125095128D+0/,W40(13)/0.648040134566009D-1/
      DATA X40(14)/0.483075801686179D+0/,W40(14)/0.679120458152339D-1/
      DATA X40(15)/0.413779204371605D+0/,W40(15)/0.706116473912869D-1/
      DATA X40(16)/0.341994090825758D+0/,W40(16)/0.728865823958039D-1/
      DATA X40(17)/0.268152185007254D+0/,W40(17)/0.747231690579679D-1/
      DATA X40(18)/0.192697580701371D+0/,W40(18)/0.761103619006259D-1/
      DATA X40(19)/0.116084070675255D+0/,W40(19)/0.770398181642479D-1/
      DATA X40(20)/0.387724175060510D-1/,W40(20)/0.775059479784249D-1/
      DATA X48( 1)/0.998771007252426D+0/,W48( 1)/0.315334605230600D-2/
      DATA X48( 2)/0.993530172266351D+0/,W48( 2)/0.732755390127600D-2/
      DATA X48( 3)/0.984124583722827D+0/,W48( 3)/0.114772345792340D-1/
      DATA X48( 4)/0.970591592546247D+0/,W48( 4)/0.155793157229440D-1/
      DATA X48( 5)/0.952987703160431D+0/,W48( 5)/0.196161604573550D-1/
      DATA X48( 6)/0.931386690706554D+0/,W48( 6)/0.235707608393240D-1/
      DATA X48( 7)/0.905879136715570D+0/,W48( 7)/0.274265097083570D-1/
      DATA X48( 8)/0.876572020274248D+0/,W48( 8)/0.311672278327980D-1/
      DATA X48( 9)/0.843588261624393D+0/,W48( 9)/0.347772225647700D-1/
      DATA X48(10)/0.807066204029443D+0/,W48(10)/0.382413510658310D-1/
      DATA X48(11)/0.767159032515740D+0/,W48(11)/0.415450829434650D-1/
      DATA X48(12)/0.724034130923815D+0/,W48(12)/0.446745608566940D-1/
      DATA X48(13)/0.677872379632664D+0/,W48(13)/0.476166584924900D-1/
      DATA X48(14)/0.628867396776514D+0/,W48(14)/0.503590355538540D-1/
      DATA X48(15)/0.577224726083973D+0/,W48(15)/0.528901894851940D-1/
      DATA X48(16)/0.523160974722233D+0/,W48(16)/0.551995036999840D-1/
      DATA X48(17)/0.466902904750958D+0/,W48(17)/0.572772921004030D-1/
      DATA X48(18)/0.408686481990717D+0/,W48(18)/0.591148396983960D-1/
      DATA X48(19)/0.348755886292161D+0/,W48(19)/0.607044391658940D-1/
      DATA X48(20)/0.287362487355455D+0/,W48(20)/0.620394231598930D-1/
      DATA X48(21)/0.224763790394689D+0/,W48(21)/0.631141922862539D-1/
      DATA X48(22)/0.161222356068892D+0/,W48(22)/0.639242385846479D-1/
      DATA X48(23)/0.970046992094629D-1/,W48(23)/0.644661644359499D-1/
      DATA X48(24)/0.323801709628690D-1/,W48(24)/0.647376968126839D-1/
      DATA X64( 1)/0.999305041735772D+0/,W64( 1)/0.178328072169600D-2/
      DATA X64( 2)/0.996340116771955D+0/,W64( 2)/0.414703326056200D-2/
      DATA X64( 3)/0.991013371476744D+0/,W64( 3)/0.650445796897800D-2/
      DATA X64( 4)/0.983336253884626D+0/,W64( 4)/0.884675982636400D-2/
      DATA X64( 5)/0.973326827789911D+0/,W64( 5)/0.111681394601310D-1/
      DATA X64( 6)/0.961008799652054D+0/,W64( 6)/0.134630478967190D-1/
      DATA X64( 7)/0.946411374858403D+0/,W64( 7)/0.157260304760250D-1/
      DATA X64( 8)/0.929569172131939D+0/,W64( 8)/0.179517157756970D-1/
      DATA X64( 9)/0.910522137078503D+0/,W64( 9)/0.201348231535300D-1/
      DATA X64(10)/0.889315445995114D+0/,W64(10)/0.222701738083830D-1/
      DATA X64(11)/0.865999398154093D+0/,W64(11)/0.243527025687110D-1/
      DATA X64(12)/0.840629296252580D+0/,W64(12)/0.263774697150550D-1/
      DATA X64(13)/0.813265315122797D+0/,W64(13)/0.283396726142590D-1/
      DATA X64(14)/0.783972358943341D+0/,W64(14)/0.302346570724020D-1/
      DATA X64(15)/0.752819907260532D+0/,W64(15)/0.320579283548510D-1/
      DATA X64(16)/0.719881850171611D+0/,W64(16)/0.338051618371420D-1/
      DATA X64(17)/0.685236313054233D+0/,W64(17)/0.354722132568820D-1/
      DATA X64(18)/0.648965471254657D+0/,W64(18)/0.370551285402400D-1/
      DATA X64(19)/0.611155355172393D+0/,W64(19)/0.385501531786160D-1/
      DATA X64(20)/0.571895646202634D+0/,W64(20)/0.399537411327200D-1/
      DATA X64(21)/0.531279464019894D+0/,W64(21)/0.412625632426230D-1/
      DATA X64(22)/0.489403145707053D+0/,W64(22)/0.424735151236530D-1/
      DATA X64(23)/0.446366017253464D+0/,W64(23)/0.435837245293230D-1/
      DATA X64(24)/0.402270157963992D+0/,W64(24)/0.445905581637560D-1/
      DATA X64(25)/0.357220158337668D+0/,W64(25)/0.454916279274180D-1/
      DATA X64(26)/0.311322871990211D+0/,W64(26)/0.462847965813140D-1/
      DATA X64(27)/0.264687162208767D+0/,W64(27)/0.469681828162100D-1/
      DATA X64(28)/0.217423643740007D+0/,W64(28)/0.475401657148300D-1/
      DATA X64(29)/0.169644420423993D+0/,W64(29)/0.479993885964580D-1/
      DATA X64(30)/0.121462819296120D+0/,W64(30)/0.483447622348030D-1/
      DATA X64(31)/0.729931217877989D-1/,W64(31)/0.485754674415030D-1/
      DATA X64(32)/0.243502926634240D-1/,W64(32)/0.486909570091400D-1/
      DATA X96( 1)/0.999689503883230D+0/,W96( 1)/0.796792065552010D-3/
      DATA X96( 2)/0.998364375863181D+0/,W96( 2)/0.185396078894692D-2/
      DATA X96( 3)/0.995981842987209D+0/,W96( 3)/0.291073181793495D-2/
      DATA X96( 4)/0.992543900323762D+0/,W96( 4)/0.396455433844469D-2/
      DATA X96( 5)/0.988054126329623D+0/,W96( 5)/0.501420274292752D-2/
      DATA X96( 6)/0.982517263563014D+0/,W96( 6)/0.605854550423596D-2/
      DATA X96( 7)/0.975939174585136D+0/,W96( 7)/0.709647079115386D-2/
      DATA X96( 8)/0.968326828463264D+0/,W96( 8)/0.812687692569876D-2/
      DATA X96( 9)/0.959688291448742D+0/,W96( 9)/0.914867123078339D-2/
      DATA X96(10)/0.950032717784437D+0/,W96(10)/0.101607705350080D-1/
      DATA X96(11)/0.939370339752755D+0/,W96(11)/0.111621020998380D-1/
      DATA X96(12)/0.927712456722308D+0/,W96(12)/0.121516046710880D-1/
      DATA X96(13)/0.915071423120898D+0/,W96(13)/0.131282295669610D-1/
      DATA X96(14)/0.901460635315852D+0/,W96(14)/0.140909417723140D-1/
      DATA X96(15)/0.886894517402420D+0/,W96(15)/0.150387210269940D-1/
      DATA X96(16)/0.871388505909296D+0/,W96(16)/0.159705629025620D-1/
      DATA X96(17)/0.854959033434601D+0/,W96(17)/0.168854798642450D-1/
      DATA X96(18)/0.837623511228187D+0/,W96(18)/0.177825023160450D-1/
      DATA X96(19)/0.819400310737931D+0/,W96(19)/0.186606796274110D-1/
      DATA X96(20)/0.800308744139140D+0/,W96(20)/0.195190811401450D-1/
      DATA X96(21)/0.780369043867433D+0/,W96(21)/0.203567971543330D-1/
      DATA X96(22)/0.759602341176647D+0/,W96(22)/0.211729398921910D-1/
      DATA X96(23)/0.738030643744400D+0/,W96(23)/0.219666444387440D-1/
      DATA X96(24)/0.715676812348967D+0/,W96(24)/0.227370696583290D-1/
      DATA X96(25)/0.692564536642171D+0/,W96(25)/0.234833990859260D-1/
      DATA X96(26)/0.668718310043916D+0/,W96(26)/0.242048417923640D-1/
      DATA X96(27)/0.644163403784967D+0/,W96(27)/0.249006332224830D-1/
      DATA X96(28)/0.618925840125468D+0/,W96(28)/0.255700360053490D-1/
      DATA X96(29)/0.593032364777572D+0/,W96(29)/0.262123407356720D-1/
      DATA X96(30)/0.566510418561397D+0/,W96(30)/0.268268667255910D-1/
      DATA X96(31)/0.539388108324357D+0/,W96(31)/0.274129627260290D-1/
      DATA X96(32)/0.511694177154667D+0/,W96(32)/0.279700076168480D-1/
      DATA X96(33)/0.483457973920596D+0/,W96(33)/0.284974110650850D-1/
      DATA X96(34)/0.454709422167743D+0/,W96(34)/0.289946141505550D-1/
      DATA X96(35)/0.425478988407300D+0/,W96(35)/0.294610899581670D-1/
      DATA X96(36)/0.395797649828908D+0/,W96(36)/0.298963441363280D-1/
      DATA X96(37)/0.365696861472313D+0/,W96(37)/0.302999154208270D-1/
      DATA X96(38)/0.335208522892625D+0/,W96(38)/0.306713761236690D-1/
      DATA X96(39)/0.304364944354496D+0/,W96(39)/0.310103325863130D-1/
      DATA X96(40)/0.273198812591049D+0/,W96(40)/0.313164255968610D-1/
      DATA X96(41)/0.241743156163840D+0/,W96(41)/0.315893307707270D-1/
      DATA X96(42)/0.210031310460567D+0/,W96(42)/0.318287588944110D-1/
      DATA X96(43)/0.178096882367618D+0/,W96(43)/0.320344562319920D-1/
      DATA X96(44)/0.145973714654896D+0/,W96(44)/0.322062047940300D-1/
      DATA X96(45)/0.113695850110665D+0/,W96(45)/0.323438225685750D-1/
      DATA X96(46)/0.812974954644249D-1/,W96(46)/0.324471637140640D-1/
      DATA X96(47)/0.488129851360490D-1/,W96(47)/0.325161187138680D-1/
      DATA X96(48)/0.162767448496020D-1/,W96(48)/0.325506144923630D-1/
!C
!C * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
!C
      IF (NPT .LE. 0) RETURN
      IF (NPT .EQ. 1) THEN
       XPT(1)= 0.0D0
       WHT(1)= 2.0D0
       RETURN
      END IF
!C
!C LOOK IF MESH IS STORED
!C
      DO I=1,MX
       IF (IOK(I) .EQ. NPT) GOTO 10
       IF (IOK(I) .LT. NPT) MXTRY=I
      END DO  
      NPT=IOK(MXTRY)
      I=MXTRY
   10 CONTINUE
!C
!C SET UP NORMAL GAUSS-LEGENDRE MESH
!C FOR NPT ODD, THE LAST (I.E. MIDDLE) TERM IS EVALUATED TWICE.
!C
      N2=(NPT+1)/2
      IC=IST(I)-1
      IOF=NPT+1
      DO I=1,N2
       XPT(I)=     -X(IC+I)
       XPT(IOF-I)=  X(IC+I)
       WHT(I)=      W(IC+I)
       WHT(IOF-I)=  W(IC+I)
      END DO
      RETURN
      END
!C
!C ***********************************************************
!C
       SUBROUTINE HARMONICS(MXPTS,NPTS,LM,RXYZ,YLM,NYLM)
C
C SPHERICAL HARMONICS GENERATOR BASED ON RECURRENCE FORMULAE
C WRITTEN BY DIRK POREZAG
C
C MXPTS: USED TO DEFINE DIMENSION OF RXYZ, YLM
C NPTS:  ACTUAL NUMBER OF POINTS
C LM:    LARGEST ANGULAR MOMENTUM REQUIRED
C RXYZ:  COORDINATES OF POINTS (normalized)
C YLM:   SPHERICAL HARMONICS AS CALCULATED ON POINTS
C NYLM:  TOTAL NUMBER OF SPHERICAL HARMONICS PER POINT
C 
C ATTENTION: THIS ROUTINE EXPECTS NORMALIZED POINTS IN ARRAY RXYZ
C
        IMPLICIT REAL*8 (A-H,O-Z)
        DIMENSION RXYZ(3,MXPTS),YLM(MXPTS,(LM+1)**2)
        DIMENSION PREFACRC(50)
        LOGICAL FIRST
        SAVE
        DATA FIRST/.TRUE./

        IF (LM .GT. 50) THEN
         PRINT *,'HARMONICS: THIS ROUTINE IS NOT STABLE FOR L > 50'
         STOP
        END IF
        NYLM=0
        IF (LM .LT. 0) RETURN
        NYLM=(LM+1)**2
        IF (NPTS .LE. 0) RETURN
C
C SET UP 1/I ARRAY
C
        IF (FIRST) THEN
         FIRST= .FALSE.
         DO I=1,50 
          PREFACRC(I)= 1.0D0/I
         END DO
        END IF
C
C L=0
C
        DO IPTS=1,NPTS
         YLM(IPTS,1)= 1.0D0
        END DO
        IF (LM .LT. 1) GOTO 200
C
C L=1
C
        DO IPTS=1,NPTS
         YLM(IPTS,2)= RXYZ(3,IPTS)
         YLM(IPTS,3)= RXYZ(1,IPTS)
         YLM(IPTS,4)= RXYZ(2,IPTS)
        END DO
C
C HIGHER L: START WITH RECURSION FOR Y_LM (M <= L-1)
C
        DO 100 L=2,LM
         IFC2L1=2*L-1
         NOW=L**2
         IL1=(L-1)**2
         IL2=(L-2)**2
         DO M=0,L-1
          FACRC=PREFACRC(L-M)
          IFCLM1=L+M-1
          IF (M .EQ. L-1) IFCLM1=0 
          AA= FACRC*IFC2L1
          BB= FACRC*IFCLM1
          NRUN=2
          IF (M .EQ. 0) NRUN=1
          DO IRUN=1,NRUN
           NOW=NOW+1
           IL1=IL1+1
           IL2=IL2+1
           DO IPTS=1,NPTS
            YLM(IPTS,NOW)= AA*YLM(IPTS,IL1)*YLM(IPTS,2)
     &                    -BB*YLM(IPTS,IL2)
           END DO
          END DO
         END DO
C
C RECURSION FOR Y_LL
C
         LAST=L**2
         NOW=(L+1)**2
         DO IPTS=1,NPTS
          YLM(IPTS,NOW)=   IFC2L1*(YLM(IPTS,LAST  )*YLM(IPTS,3)
     &                            +YLM(IPTS,LAST-1)*YLM(IPTS,4))
          YLM(IPTS,NOW-1)= IFC2L1*(YLM(IPTS,LAST-1)*YLM(IPTS,3)
     &                            -YLM(IPTS,LAST  )*YLM(IPTS,4))
         END DO
  100   CONTINUE
C
C PREFACTORS
C
  200   FOURPI=16*ATAN(1.0D0)
        ONOFPI=1.0D0/FOURPI
        DO L=0,LM
         FACFC= (2*L+1)*ONOFPI
         NOW=L**2
         DO M=0,L
          NRUN=2
          IF (M .EQ. 0) NRUN=1
          VFAC=SQRT(NRUN*FACFC)
          DO IRUN=1,NRUN
           NOW=NOW+1
           DO IPTS=1,NPTS
            YLM(IPTS,NOW)= VFAC*YLM(IPTS,NOW)
           END DO
          END DO
          FACFC= FACFC/(MAX(L-M,1)*(L+M+1))
         END DO
        END DO
        RETURN
       END SUBROUTINE
