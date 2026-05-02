      subroutine rotcoef(iwf,ispnx,nbas,athe,pcoef)
      include "PARAMS"
      include "commons.inc"
      dimension pcoef(ndh),vin(3),vout(3)

      pcoef(:)=0.0d0
      pi=4.0d0*atan(1.0d0)
      rang1=athe*pi/180.0d0
      cost=cos(rang1)
      sint=sin(rang1)

      icount=0
      do j=1,7
      do iat=1,2

        icount=icount+5

        ibas=(j-1)*10+(iat-1)*5
        PCOEF(1+ibas)=PSI_COEF(1+ibas,iwf,1,ispnx)
        PCOEF(2+ibas)=PSI_COEF(2+ibas,iwf,1,ispnx)

        px=PSI_COEF(3+ibas,iwf,1,ispnx)
        py=PSI_COEF(4+ibas,iwf,1,ispnx)
        pz=PSI_COEF(5+ibas,iwf,1,ispnx)
        
        PCOEF(3+ibas) = px
        PCOEF(4+ibas) = py*cost-pz*sint
        PCOEF(5+ibas) = py*sint+pz*cost

      enddo
      enddo

      if(icount.ne.nbas)then
              print *,'icount problem ne to nbas',icount
         call stopit
      endif
      return
      end
