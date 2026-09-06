      Subroutine GCENTER(RMSH,ngrid,ccenter)
      implicit real*8 (a-h,o-z)
      DIMENSION RMSH(3,ngrid),RMSH2(3,ngrid)
      DIMENSION tmat(3,3),ccenter(3),atomvec(3),vec(3),natomkind(20)
      CHARACTER(LEN=100) :: line
      LOGICAL :: is_element
      ALLOCATABLE :: atom1(:,:)
  
      open(100,file='lattice.data',form='formatted',status='unknown')
      rewind(100)

      read(100,'(A)') line

! Initialize number of atoms
      num_atoms = 0
! Initialize position tracker
      pos = 1

      do i = 1, len_trim(line)
! Check if the current character is a letter
         is_element = (line(i:i) >= 'A' .and. line(i:i) <= 'Z')
        ! If the current character is an element, count it as an atom
         if (is_element) then
            num_atoms = num_atoms + 1
            ! Increment the position by 6 to maintain the distance
            pos = pos + 6
         endif
        ! Update position tracker
          pos = pos + 1
      end do
!      print *, "Number of atoms in the first line: ", num_atoms

      read(100,*)
      do i=1,3
        read(100,*)(tmat(i,j),j=1,3)
      enddo

      read(100,*)
      read(100,*)(natomkind(ii),ii=1,num_atoms)
      read(100,*)

      ntotatoms=0
      do jj=1,num_atoms
        ntotatoms=ntotatoms+natomkind(jj)
      enddo

      ALLOCATE(atom1(ntotatoms,3))

      do iatom=1,ntotatoms
          read(100,*)(atom1(iatom,j), j = 1, 3)
      enddo

      ! Calculate center coordinates
      ccenter(1) = sum(atom1(:,1)) / real(ntotatoms)
      ccenter(2) = sum(atom1(:,2)) / real(ntotatoms)
      ccenter(3) = sum(atom1(:,3)) / real(ntotatoms)
      
!      print *, ccenter(1),ccenter(2),ccenter(3)
!      ! Move atoms to the origin
!      do i = 1, ntotatoms
!        atom1(i,1) = atom1(i, 1) - ccenter(1)
!        atom1(i,2) = atom1(i, 2) - ccenter(2)
!        atom1(i,3) = atom1(i, 3) - ccenter(3)
!      end do
      
!      ! Output the new atomic positions
!      print *, "New atomic positions after moving to the origin:"
!
!      do i = 1, ntotatoms
!      print *, i, ": (x, y, z) = ",atom1(i, 1), atom1(i, 2), atom1(i, 3)

!      print *,ngrid
!      do ii=1,ngrid
!        WRITE(101,*)(RMSH(jj,ii),jj=1,3)
!      enddo

       RMSH2(:,:)=0.0d0

       do igrid=1,ngrid
         do jj=1,3
!           RMSH2(jj,igrid)=RMSH(jj,igrid)-ccenter(jj)
           RMSH(jj,igrid)=RMSH(jj,igrid)+ccenter(jj)
         enddo
       enddo

      END Subroutine GCENTER
