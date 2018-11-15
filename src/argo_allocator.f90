!> @file
!! @brief This file provides a FORTRAN interface to ArgoDSM
!! @copyright Eta Scale AB. Licensed under the Eta Scale Open Source License. See the LICENSE file for details.
module argo_allocator
  use kind
  use iso_c_binding
  !
  interface
    type(c_ptr) function dynamic_alloc(s)
      use iso_c_binding
      integer(kind=8) :: s
    end function
    !
    type(c_ptr) function collective_alloc(s)
      use iso_c_binding
      integer(kind=8) :: s
    end function
  end interface
  !
contains
  subroutine alloc_4d_real(arr, s)
    integer(kind=ni) :: s(4_ni)
    real(kind=nr), pointer :: arr(:,:,:,:)
    type(c_ptr) :: cptr
    ! -----------------------------------
    !
    cptr = dynamic_alloc(8_8*product(s))
    call c_f_pointer(cptr, arr, s)
  end subroutine
  !
  subroutine coalloc_4d_real(arr, s)
    integer(kind=ni) :: s(4_ni)
    real(kind=nr), pointer :: arr(:,:,:,:)
    type(c_ptr) :: cptr
    ! -----------------------------------
    !
    cptr = collective_alloc(8_8*product(s))
    call c_f_pointer(cptr, arr, s)
  end subroutine 
  !
  subroutine free_4d_real(arr)
    real(kind=nr), pointer :: arr(:,:,:,:)
    ! -----------------------------------
    !
    call dynamic_free(arr)
  end subroutine
  !
  subroutine cofree_4d_real(arr)
    real(kind=nr), pointer :: arr(:,:,:,:)
    ! -----------------------------------
    !
    call collective_free(arr)
  end subroutine
  !
  subroutine alloc_3d_real(arr, s)
    integer(kind=ni) :: s(3_ni)
    real(kind=nr), pointer :: arr(:,:,:)
    type(c_ptr) :: cptr
    ! -----------------------------------
    !
    cptr = dynamic_alloc(8_8*product(s))
    call c_f_pointer(cptr, arr, s)
  end subroutine
  !
  subroutine coalloc_3d_real(arr, s)
    integer(kind=ni) :: s(3_ni)
    real(kind=nr), pointer :: arr(:,:,:)
    type(c_ptr) :: cptr
    ! -----------------------------------
    !
    cptr = collective_alloc(8_8*product(s))
    call c_f_pointer(cptr, arr, s)
  end subroutine
  !
  subroutine free_3d_real(arr)
    real(kind=nr), pointer :: arr(:,:,:)
    ! -----------------------------------
    !
    call dynamic_free(arr)
  end subroutine
  !
  subroutine cofree_3d_real(arr)
    real(kind=nr), pointer :: arr(:,:,:)
    ! -----------------------------------
    !
    call collective_free(arr)
  end subroutine
  !
  subroutine alloc_2d_real(arr, s)
    integer(kind=ni) :: s(2_ni)
    real(kind=nr), pointer :: arr(:,:)
    type(c_ptr) :: cptr
    ! -----------------------------------
    !
    cptr = dynamic_alloc(8_8*product(s))
    call c_f_pointer(cptr, arr, s)
  end subroutine
  !
  subroutine coalloc_2d_real(arr, s)
    integer(kind=ni) :: s(2_ni)
    real(kind=nr), pointer :: arr(:,:)
    type(c_ptr) :: cptr
    ! -----------------------------------
    !
    cptr = collective_alloc(8_8*product(s))
    call c_f_pointer(cptr, arr, s)
  end subroutine
  !
  subroutine free_2d_real(arr)
    real(kind=nr), pointer :: arr(:,:)
    ! -----------------------------------
    !
    call dynamic_free(arr)
  end subroutine
  !
  subroutine cofree_2d_real(arr)
    real(kind=nr), pointer :: arr(:,:)
    ! -----------------------------------
    !
    call collective_free(arr)
  end subroutine
end module
