! Copyright 2018 Clemens Spensberger <clemens.spensberger@uib.no>
! Permission is hereby granted, free of charge, to any person obtaining a copy of this software and associated documentation files
! (the "Software"), to deal in the Software without restriction, including without limitation the rights to use, copy, modify, merge,
! publish, distribute, sublicense, and/or sell copies of the Software, and to permit persons to whom the Software is furnished to do
! so, subject to the following conditions:
!
! The above copyright notice and this permission notice shall be included in all copies or substantial portions of the Software.
!
! THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF
! MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE
! FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN
! CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
!
module ptest
  use kind
  !
  implicit none
  !
  real(kind=nr), parameter :: nan = 0.0_nr / 0.0_nr
  logical, parameter :: grid_cyclic_ew = .true.
contains
  !
  !@ Calculates partial x derivative: ddatdx = partial(dat)/partial(x)
  !@
  !@ The routine uses 2nd-order centered differences. Returns NaN on first and last 
  !@ longitude for non-cyclic grids.
  !@
  !@ Parameters
  !@ ----------
  !@
  !@ dat : np.ndarray with shape (nz,ny,nx) and dtype float64
  !@     Data array
  !@ dx : np.ndarray with shape (ny,nx) and dtype float64
  !@     The double grid spacing in x-direction to be directly for centered differences.
  !@     ``dx(j,i)`` is expected to contain the x-distance between ``(j,i+1)`` and ``(j,i-1)``.
  !@ dy : np.ndarray with shape (ny,nx) and dtype float64
  !@     The double grid spacing in y-direction to be directly for centered differences.
  !@     ``dy(j,i)`` is expected to contain the y-distance between ``(j+1,i)`` and ``(j-1,i)``.
  !@
  !@ Other parameters
  !@ ----------------
  !@
  !@ nx : int
  !@     Grid size in x-direction.
  !@ ny : int
  !@     Grid size in y-direction.
  !@ nz : int
  !@     Grid size in z- or t-direction.
  !@
  !@ Returns
  !@ -------
  !@ np.ndarray with shape (nz,ny,nx) and dtype float64
  !@     x-derivative of ``dat``.
  !@
  !@ See Also
  !@ --------
  !@ :meth:`ddx_o4`
  subroutine ddx(res,nx,ny,nz,dat,dx,dy)
    integer(kind=ni), intent(in) :: nx,ny,nz
    real(kind=nr), intent(in)  :: dat(nz,ny,nx), dx(ny,nx), dy(ny,nx)
    real(kind=nr), intent(out) :: res(nz,ny,nx)
    !
    integer(kind=ni) :: i,j,k
    !f2py depend(nx,ny,nz) res
    !f2py depend(nx,ny) dx, dy
    ! -----------------------------------------------------------------
    !
    !$OMP PARALLEL DO
    do i = 2_ni,nx-1_ni
       do j = 1_ni,ny
          do k = 1_ni,nz
             res(k,j,i) = (dat(k,j,i+1_ni)-dat(k,j,i-1_ni))/dx(j,i)
          end do
       end do
    end do
    !$OMP END PARALLEL DO
    if (grid_cyclic_ew) then
       !$OMP PARALLEL DO
       do j = 1_ni,ny
          do k = 1_ni,nz
             res(k,j,1_ni) = (dat(k,j,  2_ni)-dat(k,j     ,nx))/dx(j,1_ni)
             res(k,j,nx  ) = (dat(k,j  ,1_ni)-dat(k,j,nx-1_ni))/dx(j,nx)
          end do
       end do
       !$OMP END PARALLEL DO
    else 
       !$OMP PARALLEL DO
       do j = 1_ni,ny
          do k = 1_ni,nz
             res(k,j,1_ni) = nan
             res(k,j,nx  ) = nan
          end do
       end do
       !$OMP END PARALLEL DO
    end if
  end subroutine
  !
  !@ Calculates partial y derivative: ddatdy = partial(dat)/partial(y)
  !@
  !@ The routine uses 2nd-order centered differences. Returns NaN on first and last 
  !@ latitude.
  !@
  !@ Parameters
  !@ ----------
  !@
  !@ dat : np.ndarray with shape (nz,ny,nx) and dtype float64
  !@     Data array
  !@ dx : np.ndarray with shape (ny,nx) and dtype float64
  !@     The double grid spacing in x-direction to be directly for centered differences.
  !@     ``dx(j,i)`` is expected to contain the x-distance between ``(j,i+1)`` and ``(j,i-1)``.
  !@ dy : np.ndarray with shape (ny,nx) and dtype float64
  !@     The double grid spacing in y-direction to be directly for centered differences.
  !@     ``dy(j,i)`` is expected to contain the y-distance between ``(j+1,i)`` and ``(j-1,i)``.
  !@
  !@ Other parameters
  !@ ----------------
  !@
  !@ nx : int
  !@     Grid size in x-direction.
  !@ ny : int
  !@     Grid size in y-direction.
  !@ nz : int
  !@     Grid size in z- or t-direction.
  !@
  !@ Returns
  !@ -------
  !@ np.ndarray with shape (nz,ny,nx) and dtype float64
  !@     y-derivative of ``dat``.
  !@
  !@ See Also
  !@ --------
  !@ :meth:`ddy_o4`
  subroutine ddy(res,nx,ny,nz,dat,dx,dy)
    integer(kind=ni), intent(in) :: nx,ny,nz
    real(kind=nr), intent(in)  :: dat(nz,ny,nx), dx(ny,nx), dy(ny,nx)
    real(kind=nr), intent(out) :: res(nz,ny,nx)
    !
    integer(kind=ni) :: i,j,k
    !f2py depend(nx,ny,nz) res
    !f2py depend(nx,ny) dx, dy
    ! -----------------------------------------------------------------
    !
    !$OMP PARALLEL DO
    do i = 1_ni,nx
       do j = 2_ni,ny-1_ni
          do k = 1_ni,nz
             res(k,j,i) = (dat(k,j+1_ni,i)-dat(k,j-1_ni,i))/dy(j,i)
          end do
       end do
    end do
    !$OMP END PARALLEL DO
    do i = 1_ni,nx
       do k = 1_ni,nz
          res(k,1_ni,i) = nan
          res(k,ny,i) = nan
       end do
    end do
  end subroutine
  !
  !@ Calculate shear deformation::
  !@
  !@     def_shear = du/dy + dv/dx
  !@
  !@ Parameters
  !@ ----------
  !@
  !@ u : np.ndarray with shape (nz,ny,nx) and dtype float64
  !@     The u-wind velocity field.
  !@ v : np.ndarray with shape (nz,ny,nx) and dtype float64
  !@     The v-wind velocity field.
  !@ dx : np.ndarray with shape (ny,nx) and dtype float64
  !@     The double grid spacing in x-direction to be directly for centered differences.
  !@     ``dx(j,i)`` is expected to contain the x-distance between ``(j,i+1)`` and ``(j,i-1)``.
  !@ dy : np.ndarray with shape (ny,nx) and dtype float64
  !@     The double grid spacing in y-direction to be directly for centered differences.
  !@     ``dy(j,i)`` is expected to contain the y-distance between ``(j+1,i)`` and ``(j-1,i)``.
  !@
  !@ Other parameters
  !@ ----------------
  !@
  !@ nx : int
  !@     Grid size in x-direction.
  !@ ny : int
  !@     Grid size in y-direction.
  !@ nz : int
  !@     Grid size in z- or t-direction.
  !@
  !@ Returns
  !@ -------
  !@ np.ndarray with shape (nz,ny,nx) and dtype float64
  !@     Calculated shear deformation.
  !@
  !@ See Also
  !@ --------
  !@ :meth:`def_stretch`, :meth:`def_total`, :meth:`def_angle`, :meth:`vor`, :meth:`div`
  subroutine def_shear(res,nx,ny,nz,u,v,dx,dy)
    integer(kind=ni), intent(in) :: nx,ny,nz
    real(kind=nr), intent(in)  :: u(nz,ny,nx), v(nz,ny,nx), dx(ny,nx), dy(ny,nx)
    real(kind=nr), intent(out) :: res(nz,ny,nx)
    !
    real(kind=nr) :: dyu(nz,ny,nx),dxv(nz,ny,nx)
    integer(kind=ni) :: i,j,k
    !f2py depend(nx,ny,nz) res, v
    !f2py depend(nx,ny) dx, dy
    ! -----------------------------------------------------------------
    !
    call ddx(dxv,nx,ny,nz,v,dx,dy)
    call ddy(dyu,nx,ny,nz,u,dx,dy)
    !
    !$OMP PARALLEL DO
    do i = 1_ni,nx
       do j = 1_ni,ny
          do k = 1_ni,nz
             res(k,j,i) = dyu(k,j,i) + dxv(k,j,i)
          end do
       end do
    end do
    !$OMP END PARALLEL DO
  end subroutine
  !
  !@ Calculate stretch deformation::
  !@
  !@    def_stretch =  du/dx - dv/dy
  !@
  !@ Parameters
  !@ ----------
  !@
  !@ u : np.ndarray with shape (nz,ny,nx) and dtype float64
  !@     The u-wind velocity field.
  !@ v : np.ndarray with shape (nz,ny,nx) and dtype float64
  !@     The v-wind velocity field.
  !@ dx : np.ndarray with shape (ny,nx) and dtype float64
  !@     The double grid spacing in x-direction to be directly for centered differences.
  !@     ``dx(j,i)`` is expected to contain the x-distance between ``(j,i+1)`` and ``(j,i-1)``.
  !@ dy : np.ndarray with shape (ny,nx) and dtype float64
  !@     The double grid spacing in y-direction to be directly for centered differences.
  !@     ``dy(j,i)`` is expected to contain the y-distance between ``(j+1,i)`` and ``(j-1,i)``.
  !@
  !@ Other parameters
  !@ ----------------
  !@
  !@ nx : int
  !@     Grid size in x-direction.
  !@ ny : int
  !@     Grid size in y-direction.
  !@ nz : int
  !@     Grid size in z- or t-direction.
  !@
  !@ Returns
  !@ -------
  !@ np.ndarray with shape (nz,ny,nx) and dtype float64
  !@     Calculated stretching deformation.
  !@
  !@ See Also
  !@ --------
  !@ :meth:`def_shear`, :meth:`def_total`, :meth:`def_angle`, :meth:`vor`, :meth:`div`
  subroutine def_stretch(res,nx,ny,nz,u,v,dx,dy)
    integer(kind=ni), intent(in) :: nx,ny,nz
    real(kind=nr), intent(in)  :: u(nz,ny,nx), v(nz,ny,nx), dx(ny,nx), dy(ny,nx)
    real(kind=nr), intent(out) :: res(nz,ny,nx)
    !
    real(kind=nr) :: dxu(nz,ny,nx),dyv(nz,ny,nx)
    integer(kind=ni) :: i,j,k
    !f2py depend(nx,ny,nz) res, v
    !f2py depend(nx,ny) dx, dy
    ! -----------------------------------------------------------------
    !
    call ddx(dxu,nx,ny,nz,u,dx,dy)
    call ddy(dyv,nx,ny,nz,v,dx,dy)
    !
    !$OMP PARALLEL DO
    do i = 1_ni,nx
       do j = 1_ni,ny
          do k = 1_ni,nz
             res(k,j,i) = dxu(k,j,i) - dyv(k,j,i)
          end do
       end do
    end do
    !$OMP END PARALLEL DO
  end subroutine
  !
  !@ Calculate the angle between the x-axis and the axis of dilatation
  !@ 
  !@ This angle goes by many different symbols:
  !@
  !@  * Spensberger and Spengler (2014): ``gamma``
  !@  * Markowski and Richardson (2011): ``alpha``
  !@  * Keyser, Reeder and Reed (1988), Lapeyre, Klein and Hua (1999): ``gamma``
  !@
  !@ Parameters
  !@ ----------
  !@
  !@ u : np.ndarray with shape (nz,ny,nx) and dtype float64
  !@     U-wind velocity.
  !@ v : np.ndarray with shape (nz,ny,nx) and dtype float64
  !@     V-wind velocity.
  !@ dx : np.ndarray with shape (ny,nx) and dtype float64
  !@     The double grid spacing in x-direction to be directly for centered differences.
  !@     ``dx(j,i)`` is expected to contain the x-distance between ``(j,i+1)`` and ``(j,i-1)``.
  !@ dy : np.ndarray with shape (ny,nx) and dtype float64
  !@     The double grid spacing in y-direction to be directly for centered differences.
  !@     ``dy(j,i)`` is expected to contain the y-distance between ``(j+1,i)`` and ``(j-1,i)``.
  !@
  !@ Other parameters
  !@ ----------------
  !@
  !@ nx : int
  !@     Grid size in x-direction.
  !@ ny : int
  !@     Grid size in y-direction.
  !@ nz : int
  !@     Grid size in z- or t-direction.
  !@
  !@ Returns
  !@ -------
  !@ np.ndarray with shape (nz,ny,nx) and dtype float64
  !@     Calculated deformation angle.
  !@
  !@ See Also
  !@ --------
  !@ :meth:`def_angle_nat`, :meth:`def_total`
  subroutine def_angle(res,nx,ny,nz,u,v,dx,dy)
    integer(kind=ni), intent(in) :: nx,ny,nz
    real(kind=nr), intent(in)  :: u(nz,ny,nx), v(nz,ny,nx), dx(ny,nx), dy(ny,nx)
    real(kind=nr), intent(out) :: res(nz,ny,nx)
    !
    real(kind=nr) :: sig_sh(nz,ny,nx), sig_st(nz,ny,nx)
    integer(kind=ni) :: i,j,k
    !f2py depend(nx,ny,nz) res, v
    !f2py depend(nx,ny) dx, dy
    ! -----------------------------------------------------------------
    !
    call def_shear(sig_sh,nx,ny,nz,u,v,dx,dy)
    call def_stretch(sig_st,nx,ny,nz,u,v,dx,dy)    
    !
    !$OMP PARALLEL DO
    do i = 1_ni,nx
       do j = 1_ni,ny
          do k = 1_ni,nz
             res(k,j,i) = 0.5_nr * atan2(sig_sh(k,j,i),sig_st(k,j,i))
          end do
       end do
    end do
    !$OMP END PARALLEL DO
  end subroutine
end module
