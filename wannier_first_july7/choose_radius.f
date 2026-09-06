!      implicit real*8 (a-h,o-z)
!      call choose_radius(16.0d0,sph)
!      print *,'sphere',sph
!      end

      subroutine choose_radius(acut,rsphere)
      implicit real*8 (a-h,o-z)
      dimension vecatm(3,100),alp(100),ncn(50),nxt(50)
      common/PTRANS/TVEC(3,3),TWIST(3,3,3),
     &       ROT(3,3),RSPHWD,NXYZ

      CHARACTER*200 line

      open(90,file="INRED")
      open(91,file="INPUT")

      do i=1,10000
         read(91,"(A)",end=10)line
         nsp=0
         do j=1,100
            if(line(j:j).ne." ")nsp=nsp+1
         enddo
         if(nsp.ne.0)then
         write(90,"(A)")trim(line)
         endif
      enddo
  10  continue

      rewind(90)
      close(91)
      
      read(90,*)nf
      
      rsphere=0.0d0

      do inf=1,nf
        read(90,*)  !ndum,ndum
        read(90,*)
        read(90,*)nat
        if(nat.gt.100)stop "choose_radius"
        do iat=1,nat
!The error is associated with the name of the SYMBOL in INPUT file (need coordinates)
         read(90,*)(vecatm(j,iat),j=1,3)
        enddo
        read(90,*)
        read(90,*)nalp
        read(90,"(A)")line
        do j=1,97
         if(line(j:j+2).eq."NUM")k=j-1
        enddo

        j=0
        do i=1,k
         if(line(i:i).ne." ")NE=i
        enddo

        j=0
        do i=k,1,-1
         if(line(i:i).ne." ")NB=i
        enddo

        lmax=0
        do i=NB,NE
         if(line(i:i).ne." ".and.line(i+1:i+1).eq." ")lmax=lmax+1
        enddo

        read(line(NB:NE),*)(ncn(l),l=1,lmax)
        read(90,*)         (nxt(l),l=1,lmax)

        do l=1,lmax
         ncn(l)=ncn(l)+nxt(l)
         if(ncn(l).gt.50)stop "choose_radius_50"
        enddo

        read(90,*)(alp(ialp),ialp=1,nalp)

        alp_min=1.0d30
        do ialp=1,nalp
         alp_min=min(alp(ialp),alp_min)
        enddo

c exp(-alp_min*rinf*rinf)=exp(-acut)
c rinf=dsqrt(acut/alp_min)

        rinf=dsqrt(acut/alp_min) !Decision

        rmax=0.0d0
         
        do iat=1,nat
         ratm=vecatm(1,iat)**2+vecatm(2,iat)**2+vecatm(3,iat)**2
         ratm=dsqrt(ratm)
         rmax=max(rmax,ratm)
         rsphere=max(rsphere,rinf+ratm)
        enddo

        print *,'rsphere and rmax:',rsphere,rmax

        do l=1,lmax
         do icon=1,ncn(l)
          read(90,*)(alp(ialp),ialp=1,nalp)
         enddo
        enddo
      enddo

      mx=0
      my=0
      mz=0

      print *,'rmax,NXYZ',rmax,NXYZ

      if(NXYZ.ge.1)mx=10
      if(NXYZ.ge.2)my=10
      if(NXYZ.ge.3)mz=10

      open(92,file="TMP")

      ntrn=0
      do nx=-mx,mx
      do ny=-my,my
      do nz=-mz,mz
       vrad=0.0d0
       do j=1,3
        vlat=nx*tvec(1,j)+ny*tvec(2,j)+nz*tvec(3,j)
        !vlat=nx*tvec(j,1)+ny*tvec(j,2)+nz*tvec(j,3)
        vrad=vrad+vlat**2
       enddo

       vrad=dsqrt(vrad)+rmax
       if(vrad.le.rsphere)then
         ntrn=ntrn+1
         write(92,*)nx,ny,nz
       endif
      enddo
      enddo
      enddo

      write(92,*)ntrn
      close(92)

      call system('tail -1 TMP >QTRANSLATIONS')
      call system('cat TMP >>QTRANSLATIONS')

      return
      end
