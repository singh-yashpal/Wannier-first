      Subroutine COULOMBPOT(RMSH,ngrid,COULOMB2)
      implicit real*8 (a-h,o-z)
      DIMENSION RMSH(3,ngrid),COULOMB2(ngrid),amesh(ngrid)
      DIMENSION atomvec(3),vec(3),atom1(3,100),zatm(100)
  

      open(99,file='XMOL.DAT',status='old')
      read(99,*)iatomtot
      read(99,*)

      atom1=0.0d0
      do i=1,iatomtot
       read(99,*)zatm(i),(atom1(ii,i),ii=1,3)
      enddo

      do igrid=1,ngrid
        coupot = 0.0d0
        amesh(:)=RMSH(:,igrid)
        do iatom=1,iatomtot
           atomvec(:)=atom1(:,iatom)
           vec(:)=atomvec(:)-amesh(:)
           denom=dsqrt(sum(vec(:)*vec(:)))
           coupot=coupot+1.0d0/denom
        enddo

        COULOMB2(igrid)=coupot 
      enddo

      RETURN
      END Subroutine COULOMBPOT

