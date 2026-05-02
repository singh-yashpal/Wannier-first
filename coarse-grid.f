      subroutine coarse_grid
      include 'PARAMS'
      include 'commons.inc'

      parameter (ifact = 3)
      parameter (igsize = ifact**3)
      parameter (alpha = 0.5d0)
!      parameter (icgrid = 20000)

      dimension sum_r(3),dummy_cm(3)
c     Coarse grid storage
!       COMMON/CGRID/RMSC(3,MAX_PTW),WMSC(MAX_PTW),NMSC
!      dimension RMSC(3,icgrid),WMSC(icgrid)

      idx_coarse = 0
      idx_fine = 0

 10   continue
      if (idx_fine + igsize - 1 .ge. NMSH) goto 20

      idx_coarse = idx_coarse + 1
      do j=1,3
         dummy_cm(j) = 0.0d0
      end do
      sum_w = 0.0d0

c     Compute center of group
      do ig = 0, igsize - 1
         icount = idx_fine + ig + 1
         do j=1,3
            dummy_cm(j)=dummy_cm(j)+RMSH(j,icount)
         end do
      end do

      do j=1,3
         dummy_cm(j)=dummy_cm(j)/dble(igsize)
      end do

c     Gaussian weighted centroid
      do j=1,3
         sum_r(j) = 0.0d0
      end do
      sum_w = 0.0d0

      do ig = 0, igsize - 1
         icount = idx_fine + ig + 1

         dx = RMSH(1, icount) - dummy_cm(1)
         dy = RMSH(2, icount) - dummy_cm(2)
         dz = RMSH(3, icount) - dummy_cm(3)
         r2 = dx*dx + dy*dy + dz*dz

         gauss_w = exp(-alpha * r2) * WMSH(icount)

         do j = 1, 3
            sum_r(j) = sum_r(j) + RMSH(j, icount) * gauss_w
         end do
         sum_w = sum_w + gauss_w
      end do

      do j = 1, 3
         RMSC(j, idx_coarse) = sum_r(j) / sum_w
      end do
      WMSC(idx_coarse) = sum_w

      idx_fine = idx_fine + igsize
      goto 10

 20   continue
      NMSC = idx_coarse

      leftover = mod(NMSH, igsize)
      if (leftover .ne. 0) then
         print *, 'Warning:', leftover,
     &           'fine grid points skipped in coarsening.'
      end if

      print *, 'coarse grid points:', NMSC

      do ii = 1, 10
         print *, (RMSC(j, ii), j = 1, 3), WMSC(ii)
      end do

      call stopit
      return
      end
