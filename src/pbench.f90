!> @file
!! @brief This file contains an ArgoDSM wrapper for module ptest
!! @copyright Eta Scale AB. Licensed under the Eta Scale Open Source License. See the LICENSE file for details.
program pbench
  use ptest
  use argo_allocator
  !
  implicit none
  !
  integer(kind=ni), parameter :: nx = 720_ni, ny = 361_ni, nz = 100_ni, nn = 20_ni
  !
  real(kind=nr), pointer :: u(:,:,:), v(:,:,:), res(:,:,:), dx(:,:), dy(:,:)
  integer(kind=ni) :: n, clck_counts_beg, clck_counts_end, clck_rate
  !--------------------------------------------------------------------
  !
  write(*,*) 'Initializing'
  call system_clock(clck_counts_beg, clck_rate)
  call argo_init(1024_8**3_8)
  !
  call coalloc_3d_real(u, [nz,ny,nx])
  call coalloc_3d_real(v, [nz,ny,nx])
  call coalloc_3d_real(res, [nz,ny,nx])
  call coalloc_2d_real(dx, [ny,nx])
  call coalloc_2d_real(dy, [ny,nx])
  !
  ! Initialisation
  dx(:,:) = 1.0_nr
  dy(:,:) = 1.0_nr
  !
  call system_clock(clck_counts_end, clck_rate)
  !
  write(*,*) 'Finished initialization, time [s]:', (clck_counts_end - clck_counts_beg) / real(clck_rate)
  write(*,*) 'Starting calculation, nn =', nn
  !
  call system_clock(clck_counts_beg, clck_rate)
  do n = 1_ni,nn
     call def_angle(res, nx,ny,nz, u,v, dx,dy)
  end do
  call system_clock(clck_counts_end, clck_rate)
  !
  write(*,*) 'Finished calculation, execution time [s]:'
  write(*,*) (clck_counts_end - clck_counts_beg) / real(clck_rate)
  !
  call cofree_3d_real(u)
  call cofree_3d_real(v)
  call cofree_3d_real(res)
  call cofree_2d_real(dx)
  call cofree_2d_real(dy)
  !
  call argo_finalize()
end program pbench
