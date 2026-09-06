
      implicit real*8(a-h,o-z)
      dimension ver(100,100)
      dimension ham(100,100),over(100,100)

      do i=1,100
      do j=1,100
      ham(i,j)=dfloat(i)*dfloat(j)
      ver(i,j)=dfloat(i)-dfloat(j)
      over(i,j)=dfloat(i)+dfloat(j)
      enddo
      enddo
      
      ver=0.0d0

      do i=1,100
      do j=1,100
         ver(j,i)=abs(ver(j,i))
         if(ver(j,i).gt.1D-10)then
         print *,'ver not zero',i,j,ver(j,i)
         endif
      enddo
      enddo


      end
