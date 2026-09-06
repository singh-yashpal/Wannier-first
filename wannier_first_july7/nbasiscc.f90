!Counts number of basis function in the central cell.
subroutine nbasiscc(ncent)
  implicit none

  integer :: ntypes
  integer :: itype,ncent
  integer :: natoms
  integer :: ns, np, nd
  integer :: ngauss_total
  integer :: ios
  character(len=256) :: line

  ngauss_total = 0

  open(unit=10, file='ISYMGEN', status='old', action='read')

  ! First line: number of inequivalent atom types
  read(10, *, iostat=ios) ntypes
  if (ios /= 0) stop 'Error reading number of atom types'

  do itype = 1, ntypes

     !--------------------------------------------------
     ! Find "NUMBER OF ATOMS OF TYPE"
     !--------------------------------------------------
     do
        read(10,'(A)',iostat=ios) line
        if (ios /= 0) stop 'EOF before NUMBER OF ATOMS OF TYPE'
        if (index(line,'NUMBER OF ATOMS OF TYPE') > 0) exit
     end do

     ! Read number of atoms of this type
     read(line,*,iostat=ios) natoms
     if (ios /= 0) stop 'Error reading number of atoms of type'

     !--------------------------------------------------
     ! Find "NUMBER OF BARE GAUSSIANS"
     !--------------------------------------------------
     do
        read(10,'(A)',iostat=ios) line
        if (ios /= 0) stop 'EOF before NUMBER OF BARE GAUSSIANS'
        if (index(line,'NUMBER OF BARE GAUSSIANS') > 0) exit
     end do

     !--------------------------------------------------
     ! Next line: primary S,P,D functions
     !--------------------------------------------------
     read(10,*,iostat=ios) ns, np, nd
     if (ios /= 0) stop 'Error reading S,P,D basis counts'

     !--------------------------------------------------
     ! Accumulate total Gaussian functions
     !--------------------------------------------------
     ngauss_total = ngauss_total + natoms * (ns + 3*np + 6*nd)

  end do

  close(10)

  ncent=ngauss_total
!  write(*,'(A,I6)') 'Total number of Gaussian basis functions = ', ngauss_total

  return
end
