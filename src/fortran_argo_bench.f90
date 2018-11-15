!> @file
!! @brief This file contains a small ArgoDSM application in FORTRAN
!! @copyright Eta Scale AB. Licensed under the Eta Scale Open Source License. See the LICENSE file for details.
program fortran_argo_bench
  use kind
  use npy
  use argo_allocator
  use ptest, only: def_angle
  !
  implicit none
  !
  real(kind=nr), pointer :: u(:,:,:,:), v(:,:,:,:), defang(:,:,:,:), &
                            &   dx(:,:), dy(:,:)
  integer(kind=ni), allocatable :: dshape(:)
  integer(kind=ni) :: uunit, vunit, dtype, tidx, &
         &            clck_counts_beg, clck_counts_end, clck_rate
  ! -------------------------------------------------------------------
  !
  call argo_init(7.5_8 * 1024_8**3_8)
  !
  call system_clock(clck_counts_beg, clck_rate)
  call open_npy(uunit, dtype, dshape, 'u.npy')
  call coalloc_4d_real(u, dshape)
  call read_f8_4d(uunit, u)
  !
  call open_npy(vunit, dtype, dshape, 'v.npy')
  call coalloc_4d_real(v, dshape)
  call read_f8_4d(vunit, v)
  !
  call close_npy(uunit)
  call close_npy(vunit)
  !
  call coalloc_4d_real(defang, dshape)
  call coalloc_2d_real(dx, dshape(3_ni:4_ni))
  call coalloc_2d_real(dy, dshape(3_ni:4_ni))
  dx(:,:) = 240.0e3
  dy(:,:) = 240.0e3
  !
  call system_clock(clck_counts_end, clck_rate)
  !
  write(*,*) 'Finished loading data, loading time [s]:'
  write(*,*) (clck_counts_end - clck_counts_beg) / real(clck_rate)
  ! The actual calculations to be benchmarked
  call system_clock(clck_counts_beg, clck_rate)
  do tidx = 1_ni,size(u,1_ni)
     call def_angle(defang(tidx,:,:,:), size(u,2_ni), size(u,3_ni), size(u,4_ni),  &
             &      u(tidx,:,:,:), v(tidx,:,:,:), dx, dy)
  end do
  call system_clock(clck_counts_end, clck_rate)
  !
  write(*,*) 'Finished calculation, execution time [s]:'
  write(*,*) (clck_counts_end - clck_counts_beg) / real(clck_rate)
  !
 end program
