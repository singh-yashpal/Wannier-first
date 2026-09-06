subroutine rewrite_symbol(ncalc,mode_run,imesh,ierr)
    implicit none
    integer, parameter :: max_atoms = 10000, max_types = 50
    character(len=100) :: line, kine, jine(100), filename
    character(len=10) :: element_type
    character(len=10), dimension(max_types) :: type_list
    real(8), dimension(3, max_atoms, max_types) :: coords, wcoords
    real(8), dimension(max_types, max_atoms) :: dist
    real(8), dimension(3) :: TV,vin,vout
    real(8), dimension(100) :: alpha, coef
    integer, dimension(3,10000) :: nt
    logical :: found, exist
    integer :: i, j, k, ii, iper, ios, icell, ic, eq_pos
    integer :: is,ip,id,iss,isp,isd,ij,ik,jk,iwx,iwy,iwz
    integer :: ishellx,ishelly,ishellz,intax,iraxis
    integer :: num_types, type_index, itotatom
    integer :: atom_counts(max_types), atom_count1(max_types)
    integer :: nxyz, nineq, nuce, nucn, nexec, iinequ
    integer :: ngauss, ncon1, ncon2, igauss
    integer :: lgmin, lgmax, lgauss,iatms
    integer :: mcalc,mode_run,imesh,ierr,ncalc,perax
    real(8) :: x, y, z, ucut, theta, RANG, pi, barcut, barsht
    real(8) :: BBOUND(3,2),TVEC(3,3)
    real(8) :: RSPHWD

!    common /PTRANS/ BBOUND, TVEC, theta, perax, RSPHWD, NXYZ

    pi=4.0d0*atan(1.0d0)

    RSPHWD=40.0d0
    tvec = 0.0d0

    open(90, file='PTRANSLATIONS', status='old', action='read', iostat=ios)
    if (ios /= 0) then
        print *, "Error opening file PTRANSLATIONS"
        call stopit
    end if

    read(90,*)NXYZ
    do i=1,3
    read(90,*)BBOUND(i,1),BBOUND(i,2)
    enddo
    do iper=1,3
    read(90,*)(TVEC(j,iper),j=1,3)
    end do
    read(90,*)ucut
    read(90,*)theta,iraxis
    read(90,*)iwx,iwy,iwz
    close(90)

! generating UTRANS using number of shells "ishell".
    open(unit=9, file='UTRANS')
    rewind(9)

    do i=-iwx,iwx
    do j=-iwy,iwy
    do k=-iwz,iwz
       if(i.eq.0 .and.j.eq.0.and.k.eq.0)goto 91
       write(9,*)i,j,k
 91 continue
    enddo
    enddo
    enddo
    close(9)

    INQUIRE(FILE='ISYMTOT1',EXIST=EXIST)
    IF (EXIST) GOTO 999

    ! Initialize
    num_types = 0
    atom_counts = 0
    atom_count1 = 0
    type_list = ""
    coords = 0.0d0
    wcoords = 0.0d0

    open(unit=10, file='SYMBOL', status='old', action='read', iostat=ios)
    if (ios /= 0) then
        print *, "Error opening file SYMBOL"
        call stopit
    end if

    ! Begin reading atom entries
    do
        read(10, '(A)', iostat=ios) line
        if (ios /= 0) exit
        if (index(line, 'ALL-') == 1) then
            ! Extract type from 'ALL-XXX###'
            element_type = line(5:7)

            ! Extract coordinates
            eq_pos = index(line, "=")
            if (eq_pos > 0) then
                read(line(eq_pos+1:), *) x, y, z
            else
                print *, "Malformed line: ", trim(line)
                cycle
            end if

            ! Check if this type already exists
            found = .false.
            do j = 1, num_types
                if (trim(type_list(j)) == trim(element_type)) then
                    type_index = j
                    found = .true.
                    exit
                end if
            end do

! Add new type if not found
            if (.not. found) then
                num_types = num_types + 1
                type_index = num_types
                type_list(type_index) = element_type
            end if

            ! Store coordinates
            atom_count1(type_index) = atom_count1(type_index) + 1
            coords(1, atom_count1(type_index), type_index) = x
            coords(2, atom_count1(type_index), type_index) = y
            coords(3, atom_count1(type_index), type_index) = z
        end if
    end do

    close(10)

    ! Print results
    print *, "Atom types and counts:"
    do i = 1, num_types
        print *, trim(type_list(i)), ": ", atom_count1(i), " atoms"
        do j = 1, atom_count1(i)
            print '(A10, 3F12.6)', trim(type_list(i)), coords(1,j,i), coords(2,j,i), coords(3,j,i)
        end do
    end do

! UTRANSManual
    INQUIRE(FILE='UTRANSM',EXIST=EXIST)
    IF (EXIST) CALL system('cp UTRANSM UTRANS')
! Now reading UTRANS and get coordinates of the neighbouring cell atoms.
    nt=0
    open(20, file='UTRANS', status='old', action='read', iostat=ios)
    if (ios /= 0) then
        print *, "Error opening file UTRANS"
        call stopit
    end if

! Try reading one line to check if file is empty
    rewind(20)
    read(20, *, iostat=ios) nt(1,1), nt(2,1), nt(3,1)
    if (ios /= 0) then
    print *, "UTRANS file is empty!"
    close(20)
    call stopit
    end if

    ! Start reading from beginning
    rewind(20)
    icell = 1
    do
        read(20, *, iostat=ios) nt(1,icell), nt(2,icell), nt(3,icell)
        if (ios /= 0) exit
        icell = icell + 1
    end do

    print *,'Total cells in Wannier Domain:',icell
    do i=1,icell-1
      print *,'cell',(nt(j,i),j=1,3)
    enddo
    close(20)

    iatms=0
    do j=1,num_types
    do i=1,atom_count1(j)
     iatms=iatms+1
    enddo
    enddo

    open(401,file='CELL.DAT',status='unknown')
    write(401,*)iatms
    write(401,*)
    do j=1,num_types
    do i=1,atom_count1(j)
     write (401,*)(coords(k,i,j),k=1,3)
    enddo
    enddo
    close(401)
  
!    print *, "Translated atom positions for neighboring cells:"
    do j = 1, num_types  ! Loop over atom types
       ic=0
       do ii = 1, atom_count1(j)  ! Loop over atoms in that type
          ! Coordinates in central cell
          x = coords(1,ii,j)
          y = coords(2,ii,j)
          z = coords(3,ii,j)
       do i = 1,icell-1  ! Loop over UTRANS entries
          ic = ic + 1
          do k=1,3
          TV(k)=nt(1,i)*TVEC(k,1)+nt(2,i)*TVEC(k,2)+nt(3,i)*TVEC(k,3)
          enddo
          ! Translation: only x-direction matters in 1D system
          vin(1) = x + TV(1)
          vin(2) = y + TV(2)
          vin(3) = z + TV(3)

          if(iraxis.eq.1)RANG=nt(1,i)*theta        
          if(iraxis.eq.2)RANG=nt(2,i)*theta        
          if(iraxis.eq.3)RANG=nt(3,i)*theta        

          if(abs(RANG).gt.1.0d-10)then
          call rotvec(iraxis,RANG,vin,vout)
          else
          vout(1)=vin(1)
          vout(2)=vin(2)
          vout(3)=vin(3)
          endif
          wcoords(1,ic,j) = vout(1)
          wcoords(2,ic,j) = vout(2)
          wcoords(3,ic,j) = vout(3)

       end do
       end do
       atom_counts(j)=ic
    end do

    !print *, "Final atom counts including neighbors:"
    itotatom=0
    do i = 1, num_types
      print '(A10,2X, I2)',trim(type_list(i)), atom_counts(i)
      itotatom=itotatom+atom_count1(i)+atom_counts(i)
    end do
    !print *,'TOTAL ATOMS IN Wannier Domain:',itotatom

    do j=1,num_types
    do i=1,atom_counts(j)
       dist(j,i)=dsqrt(wcoords(1,i,j)**2+wcoords(2,i,j)**2+wcoords(3,i,j)**2)
       print '(A10, 3F12.6)', trim(type_list(j)),(coords(k,i,j),k=1,3)
    enddo
    enddo
 
!calculating distance from origin

! Bubble sort for each atom type
do j = 1, num_types   ! Loop over atom types
    do i = 1, atom_counts(j) - 1
        do k = i + 1, atom_counts(j)
            if (dist(j,i) > dist(j,k)) then
                ! Swap distances
                call swap_real(dist(j,i), dist(j,k))

                ! Swap coordinates (x, y, z)
                call swap_vec(wcoords(:,i,j), wcoords(:,k,j))
            end if
        end do
    end do
end do

open(30,file='ALLCELLS')

!writing xyz file with atoms in the wannier domain including central cell.
open(28,file='XMOLWD.DAT')

open(27,file='XMOL1.DAT')
write(27,*)itotatom
write(27,*)'atoms in wannier domain'

open(29,file='XMOLWD.xyz')
write(29,*)itotatom
write(29,*)'atoms in wannier domain'

do j=1,num_types
    do i=1,atom_count1(j)
        write(29,'(A3,3F12.6)')trim(type_list(j)),(0.5291772*coords(k,i,j),k=1,3)
        write(28,'(A3,3F12.6)')trim(type_list(j)),(coords(k,i,j),k=1,3)
        write(27,'(A,3F12.6)')'6',(0.5291772*coords(k,i,j),k=1,3)
    end do
end do

do j=1,num_types
    do i=1,atom_counts(j)
!       !print '(A10, 3F12.6)', trim(type_list(j)),(wcoords(k,i,j),k=1,3)
        write(27,'(A,3F12.6)')'6',(0.5291772*wcoords(k,i,j),k=1,3)
        write(29,'(A3,3F12.6)')trim(type_list(j)),(0.5291772*wcoords(k,i,j),k=1,3)
        write(28,'(A3,3F12.6)')trim(type_list(j)),(wcoords(k,i,j),k=1,3)
        write(30,'(A3,I1,A3,I3.3,1X,A1,7X,3F12.6)') 'ALL', j, trim(type_list(j)), i, '=', (wcoords(k,i,j),k=1,3)
    end do
end do

close(30)
close(29)
close(28)
close(27)

!create SYMTOT filel
    print *,'Creating SYMTOT file'
    call system('rm SYMTOT')

    call system('head -n5 SYMBOL > SYMTOT')

    open(31,file='SYMHLP')
    write(31,'(I3,1X,A30)') (2+itotatom), 'NUMBER OF SYMBOLS IN LIST'
    write(31,'(I3,1X,A30)') itotatom, 'NUMBER OF NUCLEI'
    do i=1,itotatom
    write(31,*)'1.0 1 1 1'
    enddo
    write(31,*)'  1  CALCULATION SETUP BY ISETUP'
    close(31)

    call system('grep "^ALL-" SYMBOL >> SYMHLP')
    call system('cat SYMHLP >>SYMTOT')
    call system('cat ALLCELLS >>SYMTOT')
    call system('tail -n2 SYMBOL >>SYMTOT')
!    call system("if [ $(uname) = Darwin ]; then sed -i '''5s/ISYMGEN = INPUT/ISYMTOT = INPUT/' SYMTOT; else sed -i '5s/ISYMGEN = INPUT/ISYMTOT = INPUT/' SYMTOT; fi")

     call fix_symtot()

!    call system('sed -i '''' "5s/ISYMGEN = INPUT/ISYMTOT = INPUT/" SYMTOT') !For MAC
!    call system('sed -i "5s/ISYMGEN = INPUT/ISYMTOT = INPUT/" SYMTOT') !For Linux

! Creating ISYMTOT file
    call system('awk ''NR>1 { a[NR]=$0 } END { for (i=2; i<NR-1; i++) print a[i] }'' ISYMGEN > ISYMADD')
    call system('rm ISYMTOT')

!    do j=1,num_types
!       print '(A10)', trim(type_list(j))
!    enddo


    open(40,file='ISYMGEN')
    rewind(40)
    read(40,*)NINEQ

    open(39,file='ISYMTOT')
    write(39,'(4X,I1,3X,A)')NINEQ+num_types,'!NUMBER OF FUNCTION SETS'
    close(39)
    call system('cat ISYMADD>>ISYMTOT')
    call system('rm ISYMADD')

    DO 111 k=1,num_types
       rewind(40)
       OPEN(45,file='ISYMADD')
       DO i=1,10000
         READ(40,'(A)',END=15)KINE
         DO J=1,50
         IF(KINE(J:J+6).EQ.'ATOMS O')READ(KINE,*)NINEQ
         IF(KINE(J:J+3).EQ.trim(type_list(k)))THEN
               !READ(KINE,*)NINEQ
               !print *,'number of inequivalent atoms',NINEQ
               write(45,*)'   0   0      ELECTRONIC AND NUCLEAR CHARGE'
               write(45,'(A)')'ALL'
               write(45,'(I5,1X,A,I0,A,A)') atom_counts(k), 'NUMBER OF ATOMS OF TYPE:ALL', k, trim(type_list(k))
               GOTO 16
         ENDIF
         ENDDO
       ENDDO
  16   CONTINUE

       do j=1,atom_counts(k)
       write(45,'(A,I0,A,I3.3)')'ALL',k,trim(type_list(k)),j
       enddo
       write(45,'(A)')'EXTRABASIS'

       DO j=1,NINEQ+1
       READ(40,*)
       ENDDO

       READ(40,*)NGAUSS
       READ(40,*)is,ip,id
       READ(40,*)iss,isp,isd
       WRITE(45,'(I5,2X,A)')ngauss, ' NUMBER OF BARE GAUSSIANS'
       WRITE(45,'(3I5,2X,A)')is,ip,id, 'NUMBER OF S, P, D FUNCTIONS'
       WRITE(45,'(3I5,2X,A)')iss,isp,isd, 'SUPPLEMENTRY S, P, D FUNCTIONS'

       READ(40,*)(ALPHA(IGAUSS),IGAUSS=1,NGAUSS)
       WRITE(45,'(3D20.10)')(ALPHA(IGAUSS),IGAUSS=1,NGAUSS)
       WRITE(45,*)' '
       do j=1,is+iss
       READ(40,*)(COEF(igauss),igauss=1,ngauss)
       WRITE(45,'(3D20.10)')(COEF(igauss),igauss=1,ngauss)
       WRITE(45,*)' '
       enddo
       do j=1,ip+isp
       READ(40,*)(COEF(igauss),igauss=1,ngauss)
       WRITE(45,'(3D20.10)')(COEF(igauss),igauss=1,ngauss)
       WRITE(45,*)' '
       enddo
       do j=1,id+isd
       READ(40,*)(COEF(igauss),igauss=1,ngauss)
       WRITE(45,'(3D20.10)')(COEF(igauss),igauss=1,ngauss)
       WRITE(45,*)' '
       enddo

       call system('cat ISYMADD>>ISYMTOT')
       close(45)
       call system('rm ISYMADD')
       
 111   CONTINUE
  15   CONTINUE

       close(40)
       call system('tail -n2 ISYMGEN >ISYMTEMP')
       call system('cat ISYMTEMP >>ISYMTOT')

  999   CONTINUE
        print*,'END OF REWRITE_SYMBOL'
        call system('cp ISYMTOT1 ISYMTOT')
!        call system('cp INPUT1 INPUT')
        call system('cp SYMTOT1 SYMTOT')

        ierr=-1
       call symbol(mcalc,mode_run,imesh,ierr)
       !call stopit
       return
end 

!======================================================================================
    subroutine swap_real(a, b)
        real(8), intent(inout) :: a, b
        real(8) :: tmp,ux,uy,uz
        tmp = a
        a = b
        b = tmp
    end subroutine swap_real

    subroutine swap_vec(vec1, vec2)
        real(8), intent(inout) :: vec1(3), vec2(3)
        real(8) :: tmp(3)
        tmp = vec1
        vec1 = vec2
        vec2 = tmp
    end subroutine swap_vec



    subroutine rotvec(pax, ang1, vin, vout)
    implicit none

    integer, intent(in) :: pax
    real(8), intent(in) :: ang1
    real(8), intent(in) :: vin(3)
    real(8), intent(out) :: vout(3)

    real(8) :: rmat(3,3)
    real(8) :: ang, ux, uy, uz, pi
    integer :: i, j

    pi = 4.0d0 * atan(1.0d0)
    ang = ang1 * (pi / 180.0d0)

    rmat = 0.0d0
    ux = 0.0d0
    uy = 0.0d0
    uz = 0.0d0

    select case (pax)
    case (1)
        ux = 1.0d0
    case (2)
        uy = 1.0d0
    case (3)
        uz = 1.0d0
    case default
        print *, 'Error in rotvec: invalid rotation axis pax =', pax
        stop
    end select

    rmat(1,1) = cos(ang) + ux*ux*(1.0d0 - cos(ang))
    rmat(1,2) = ux*uy*(1.0d0 - cos(ang)) - uz*sin(ang)
    rmat(1,3) = ux*uz*(1.0d0 - cos(ang)) + uy*sin(ang)

    rmat(2,1) = uy*ux*(1.0d0 - cos(ang)) + uz*sin(ang)
    rmat(2,2) = cos(ang) + uy*uy*(1.0d0 - cos(ang))
    rmat(2,3) = uy*uz*(1.0d0 - cos(ang)) - ux*sin(ang)

    rmat(3,1) = uz*ux*(1.0d0 - cos(ang)) - uy*sin(ang)
    rmat(3,2) = uz*uy*(1.0d0 - cos(ang)) + ux*sin(ang)
    rmat(3,3) = cos(ang) + uz*uz*(1.0d0 - cos(ang))

    do i = 1, 3
        vout(i) = 0.0d0
        do j = 1, 3
            vout(i) = vout(i) + rmat(i,j) * vin(j)
        end do
    end do

    end subroutine rotvec

!    subroutine rotvec(pax,ang1,vin,vout)
!        real(8), dimension(3,3) :: rmat
!        real(8), dimension(3) :: vin,vout
!        integer :: pax,i,j
!        real(8) :: ang,ang1,ux,uy,uz,pi
!
!         pi=4.0d0*atan(1.0d0)
!         RMAT=0.0d0
!         if(pax.eq.1)then
!         ux=1.0d0
!         uy=0.0d0
!         uz=0.0d0
!         endif
!         if(pax.eq.2)then
!         ux=0.0d0
!         uy=1.0d0
!         uz=0.0d0
!         endif
!         if(pax.eq.3)then
!         ux=0.0d0
!         uy=0.0d0
!         uz=1.0d0
!         endif
!
!         ang=ang1*(pi/180.0d0)
!         RMAT(1,1)=cos(ang)+ux**2*(1.0d0-cos(ang))
!         RMAT(1,2)=ux*uy*(1.0d0-cos(ang))-uz*sin(ang)
!         RMAT(1,3)=ux*uz*(1.0d0-cos(ang))+uy*sin(ang)
!                  
!         RMAT(2,1)=uy*ux*(1.0d0-cos(ang))+uz*sin(ang)
!         RMAT(2,2)=cos(ang)+uy**2*(1.0d0-cos(ang))
!         RMAT(2,3)=uy*uz*(1.0d0-cos(ang))-ux*sin(ang)
!                  
!         RMAT(3,1)=ux*uz*(1.0d0-cos(ang))-uy*sin(ang)
!         RMAT(3,2)=uy*uz*(1.0d0-cos(ang))+ux*sin(ang)
!         RMAT(3,3)=cos(ang)+uz**2*(1.0d0-cos(ang))
!                 
!         !print *,'ang,ang1',ang,ang1
!!         do i=1,3
!!            print *,(RMAT(j,i),j=1,3)
!!         enddo
!
!         do i = 1, 3
!            vout(i) = 0.0d0
!            do j = 1, 3
!               vout(i) = vout(i) + RMAT(i,j) * vin(j)
!            end do
!         end do
!         return
!    end subroutine rotvec

    subroutine fix_symtot()
    implicit none
    character(len=256) :: cmd

  ! Construct shell command to replace line 5 in SYMTOT depending on OS
    cmd = "if [ $(uname) = Darwin ]; then " // &
        "sed -i '' '5s/ISYMGEN = INPUT/ISYMTOT = INPUT/' SYMTOT; " // &
        "else sed -i '5s/ISYMGEN = INPUT/ISYMTOT = INPUT/' SYMTOT; fi"

    call system(trim(cmd))
    end subroutine fix_symtot
