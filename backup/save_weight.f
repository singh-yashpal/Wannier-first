      subroutine save_rhow(NMSHW)
      include 'PARAMS'
      include 'commons.inc'

      character*20 fname
      integer iunit, iiw
      save iiw
      data iiw /1/

      iunit = 402

c     Create filename like 'RHOW_1', 'RHOW_2', ...
      write(fname, '("RHOW_",I3.3)') iiw

      open(iunit, file=fname, status='unknown', form='formatted')
      rewind(iunit)

      NGRAD=1
      do IMSH = 1, NMSHW
        do ISPN = 1, NSPN
          do IGRAD = 1, NGRAD
            write(iunit, *) WMSH(IMSH)
c           Optionally: write(iunit, *) RHOG(IMSH,IGRAD,ISPN), WMSH(IMSH)
          end do
        end do
      end do

      close(iunit)

c     Increment for next call
      iiw = iiw + 1

      return
      end

