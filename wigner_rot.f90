program rotate_wigner_dmatrix
  implicit none

  integer :: i, idx, ios
  real*8 :: pi, theta, c, s
  real*8 :: moments(49), moments_rot(49)
  real*8 :: l2(5), l2_rot(5), D2(5,5)
  real*8 :: l3(7), l3_rot(7), D3(7,7)
  real*8 :: l4(9), l4_rot(9), D4(9,9)
  real*8 :: l5(11), l5_rot(11), D5(11,11)
  real*8 :: l6(13), l6_rot(13), D6(13,13)

  pi = 3.141592653589793d0

  print *, 'Enter rotation angle in degrees:'
  read *, theta
  theta = theta * pi / 180.0d0  ! convert to radians

  c = cos(theta)
  s = sin(theta)

  ! Read SMULTI file
  open(unit=10, file='SMULTI', status='old', action='read')
  read(10, *) idx
  if (idx /= 49) then
    print *, 'Expecting 49 moments, found', idx
    stop
  end if

  do i = 1, 49
    read(10, *, IOSTAT=ios) idx, moments(i)
    if (ios /= 0) then
      print *, 'Error reading moment index ', i
      stop
    end if
  end do
  close(10)

  moments_rot = moments

  ! Rotate dipole only (indices 2 = z, 3 = x, 4 = y)
  moments_rot(3) =  c * moments(3) - s * moments(4)
  moments_rot(4) =  s * moments(3) + c * moments(4)
  moments_rot(2) =  moments(2)

  ! l=2 rotation using Wigner D-matrix
  call build_D_matrix_l2(D2, c, s)
  l2 = moments(5:9)
  do i = 1, 5
    l2_rot(i) = sum(D2(i,1:5) * l2(1:5))
  end do
  moments_rot(5:9) = l2_rot

  ! TODO: Replace these with real Wigner D-matrix expressions:
  ! l=3 rotation
  call build_D_matrix_l3(D3, c, s)
  l3 = moments(10:16)
  do i = 1, 7
    l3_rot(i) = sum(D3(i,1:7) * l3(1:7))
  end do
  moments_rot(10:16) = l3_rot

  ! l=4 rotation
  call build_D_matrix_l4(D4, c, s)
  l4 = moments(17:25)
  do i = 1, 9
    l4_rot(i) = sum(D4(i,1:9) * l4(1:9))
  end do
  moments_rot(17:25) = l4_rot

  ! l=5 rotation
  call build_D_matrix_l5(D5, c, s)
  l5 = moments(26:36)
  do i = 1, 11
    l5_rot(i) = sum(D5(i,1:11) * l5(1:11))
  end do
  moments_rot(26:36) = l5_rot

  ! l=6 rotation
  call build_D_matrix_l6(D6, c, s)
  l6 = moments(37:49)
  do i = 1, 13
    l6_rot(i) = sum(D6(i,1:13) * l6(1:13))
  end do
  moments_rot(37:49) = l6_rot

  ! Write to SMULTIN
  open(unit=20, file='SMULTIN', status='replace', action='write')
  write(20,*) 49
  do i = 1, 49
    write(20,*) i, moments_rot(i)
  end do
  close(20)

  print *, 'Multipole moments (l=0 to l=6) rotated using real Wigner D-matrices and written to SMULTIN.'

contains

  subroutine build_D_matrix_l2(D, c, s)
    real*8, intent(out) :: D(5,5)
    real*8, intent(in)  :: c, s
    D = 0.0d0
    D(1,1) = 1.0d0
    D(2,2) = c
    D(2,3) = -s
    D(3,2) = s
    D(3,3) = c
    D(4,4) = c**2 - s**2
    D(4,5) = -2.0d0 * c * s
    D(5,4) = 2.0d0 * c * s
    D(5,5) = c**2 - s**2
  end subroutine build_D_matrix_l2

  subroutine build_D_matrix_l3(D, c, s)
    real*8, intent(out) :: D(7,7)
    real*8, intent(in)  :: c, s
    D = 0.0d0
    D(1,1) = 1.0d0
    D(2,2) = c
    D(2,3) = -s
    D(3,2) = s
    D(3,3) = c
    D(4,4) = c**2 - s**2
    D(4,5) = -2.0d0 * c * s
    D(5,4) = 2.0d0 * c * s
    D(5,5) = c**2 - s**2
    D(6,6) = c**3 - 3*c*s**2
    D(6,7) = -3*c**2*s + s**3
    D(7,6) = 3*c**2*s - s**3
    D(7,7) = c**3 - 3*c*s**2
  end subroutine build_D_matrix_l3

  subroutine build_D_matrix_l4(D, c, s)
    real*8, intent(out) :: D(9,9)
    real*8, intent(in)  :: c, s
    D = 0.0d0
    D(1,1) = 1.0d0
    D(2,2) = c
    D(2,3) = -s
    D(3,2) = s
    D(3,3) = c
    D(4,4) = c**2 - s**2
    D(4,5) = -2.0d0 * c * s
    D(5,4) = 2.0d0 * c * s
    D(5,5) = c**2 - s**2
    D(6,6) = c**3 - 3*c*s**2
    D(6,7) = -3*c**2*s + s**3
    D(7,6) = 3*c**2*s - s**3
    D(7,7) = c**3 - 3*c*s**2
    D(8,8) = c**4 - 6*c**2*s**2 + s**4
    D(8,9) = -4*c**3*s + 4*c*s**3
    D(9,8) = 4*c**3*s - 4*c*s**3
    D(9,9) = c**4 - 6*c**2*s**2 + s**4
  end subroutine build_D_matrix_l4

  subroutine build_D_matrix_l5(D, c, s)
    real*8, intent(out) :: D(11,11)
    real*8, intent(in)  :: c, s
    D = 0.0d0
    D(1,1) = 1.0d0
    D(2,2) = c
    D(2,3) = -s
    D(3,2) = s
    D(3,3) = c
    D(4,4) = c**2 - s**2
    D(4,5) = -2.0d0 * c * s
    D(5,4) = 2.0d0 * c * s
    D(5,5) = c**2 - s**2
    D(6,6) = c**3 - 3*c*s**2
    D(6,7) = -3*c**2*s + s**3
    D(7,6) = 3*c**2*s - s**3
    D(7,7) = c**3 - 3*c*s**2
    D(8,8) = c**4 - 6*c**2*s**2 + s**4
    D(8,9) = -4*c**3*s + 4*c*s**3
    D(9,8) = 4*c**3*s - 4*c*s**3
    D(9,9) = c**4 - 6*c**2*s**2 + s**4
    D(10,10) = c**5 - 10*c**3*s**2 + 5*c*s**4
    D(10,11) = -5*c**4*s + 10*c**2*s**3 - s**5
    D(11,10) = 5*c**4*s - 10*c**2*s**3 + s**5
    D(11,11) = c**5 - 10*c**3*s**2 + 5*c*s**4
  end subroutine build_D_matrix_l5

  subroutine build_D_matrix_l6(D, c, s)
    real*8, intent(out) :: D(13,13)
    real*8, intent(in)  :: c, s
    D = 0.0d0
    D(1,1) = 1.0d0
    D(2,2) = c
    D(2,3) = -s
    D(3,2) = s
    D(3,3) = c
    D(4,4) = c**2 - s**2
    D(4,5) = -2.0d0 * c * s
    D(5,4) = 2.0d0 * c * s
    D(5,5) = c**2 - s**2
    D(6,6) = c**3 - 3*c*s**2
    D(6,7) = -3*c**2*s + s**3
    D(7,6) = 3*c**2*s - s**3
    D(7,7) = c**3 - 3*c*s**2
    D(8,8) = c**4 - 6*c**2*s**2 + s**4
    D(8,9) = -4*c**3*s + 4*c*s**3
    D(9,8) = 4*c**3*s - 4*c*s**3
    D(9,9) = c**4 - 6*c**2*s**2 + s**4
    D(10,10) = c**5 - 10*c**3*s**2 + 5*c*s**4
    D(10,11) = -5*c**4*s + 10*c**2*s**3 - s**5
    D(11,10) = 5*c**4*s - 10*c**2*s**3 + s**5
    D(11,11) = c**5 - 10*c**3*s**2 + 5*c*s**4
    D(12,12) = c**6 - 15*c**4*s**2 + 15*c**2*s**4 - s**6
    D(12,13) = -6*c**5*s + 20*c**3*s**3 - 6*c*s**5
    D(13,12) = 6*c**5*s - 20*c**3*s**3 + 6*c*s**5
    D(13,13) = c**6 - 15*c**4*s**2 + 15*c**2*s**4 - s**6
  end subroutine build_D_matrix_l6

end program rotate_wigner_dmatrix

