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
! numpy file parser for npy versions 1.0 and 2.0,
! but only supporting fortran-contiguous int16, int32 or float64 arrays
module npy
  use kind
  !
  implicit none
  !
  integer(kind=nb), parameter :: magic_npy(6_ni) = (/ -109_nb, &
     & ichar('N', kind=nb), ichar('U', kind=nb), ichar('M', kind=nb), ichar('P', kind=nb), ichar('Y', kind=nb) /)
  integer(kind=ni), parameter :: runit0 = 23_ni, max_files=20_ni
  !
  integer(kind=ni) :: startpos(max_files) = -1_ni
  integer(kind=ni) :: next_runit = runit0
contains
  subroutine open_npy(runit, dtype, dshape, filename)
    character(len=*), intent(in) :: filename
    integer(kind=ni), allocatable, intent(out) :: dshape(:)
    integer(kind=ni), intent(out) :: runit, dtype
    !
    integer(kind=nb) :: magic(6_ni), ver_major, ver_minor, n, dim_
    integer(kind=ns) :: header_len2
    integer(kind=ni) :: header_len, meta_pos, pa,pz
    character(len=:), allocatable :: header, dict, key, val
    ! -----------------------------------------------------------------
    !
    if ( next_runit >= runit0 + max_files ) then
       write(*,*) 'Error reading file: Maximum number of files (', max_files, ') reached'
       return
    end if
    !
    runit = next_runit
    next_runit = next_runit + 1_ni
    dtype = -1_ni
    !
    open(unit=runit, file=filename, status='old', access='stream')
    !
    ! Read magic bytes to identify file type
    read(runit, pos=1_ni) magic
    if ( any(magic /= magic_npy) ) then
       write(*,*) 'Error reading file: Filetype not recognized'
       write(*,*) 'File header: ', magic
       write(*,*) 'Expected header: ', magic_npy
       return
    end if
    !
    ! Read file version
    read(runit, pos=7_ni) ver_major
    read(runit, pos=8_ni) ver_minor
    if ( ver_major < 1_nb .or. ver_major > 2_nb .or. ver_minor /= 0_nb ) then
       write(*,*) 'Error reading file: Unrecognized file type version'
       write(*,*) 'Found version: ', ver_major, '.', ver_minor
       return
    end if
    ! 
    ! Read header length
    if ( ver_major == 2_nb ) then
       ! 4-bytes header length
       read(runit, pos=9_ni) header_len
       meta_pos = 13_ni
    else 
       ! 2-bytes header length, cast to 4-bytes
       read(runit, pos=9_ni) header_len2
       header_len = header_len2
       meta_pos = 11_ni
    end if
    !
    ! Read metadata
    allocate(character(header_len) :: header, dict, key, val)
    read(runit, pos=meta_pos) header
    !
    ! Set start position of the actual data
    startpos(runit-runit0+1_ni) = header_len + meta_pos
    !
    ! Parse the metadata (pretty-printed python dictionary)
    ! Note: Won't work with nested dictionaries!
    pa = index(header, '{')
    pz = index(header, '}', .true.)
    dict = header(pa+1_ni:pz-1_ni)
    !
    ! Iterating backwards over key-value pairs
    ! Backwards, because kommas might appear as part of val, but not as part of the key.
    ! Hence indexing locations (colons before vals, commas before keys) are uniquely 
    ! defined without needing to parse the substructure of val
    pz = index(dict, ',', .true.)
    dict = dict(1_ni:pz-1_ni)
    do while ( pz > 0_ni )
       pa = index(dict, ':', .true.)
       val = dict(pa+2_ni:pz-1_ni)
       !
       dict(pa:) = ' '
       !
       pz = index(dict, ',', .true.)
       if ( pz < 1_ni ) pz = -1_ni
       key = dict(pz+2_ni:pa-1_ni)
       !
       if ( pz > 0_ni ) then
          dict(pz:) = ' '
       end if
       !
       select case ( trim(key) ) 
          case ( "'descr'" )
             if ( val(2_ni:2_ni) /= '<' ) then
                write(*,*) 'Error reading file: Data must be stored in little-endian'
                return
             end if
             if ( trim(val(3_ni:4_ni)) == 'f8' ) then
                dtype = nr
             else if ( trim(val(3_ni:4_ni)) == 'i4' ) then
                dtype = ni
             else if ( trim(val(3_ni:4_ni)) == 'i2') then
                dtype = ns
             else
                write(*,*) 'Error reading file: Unknown dtype', val(3_ni:4_ni)
                return
             end if
          case ( "'fortran_order'" )
             if ( trim(val) /= 'True' ) then
                write(*,*) 'Error reading file: Data must be Fortran-contiguous'
                return
             end if
          case ( "'shape'" )
             ! Count number of commas to be able to allocate dshape
             dim_ = 1_nb
             do n = 1_nb,len(val, kind=nb)
                if ( val(n:n) == ',' ) dim_ = dim_ + 1_nb
             end do
             allocate(dshape(dim_))
             !
             val(1_ni:1_ni) = ' '
             val(len_trim(val):) = ' '
             pa = index(val, ',')
             n = 1_nb
             do while ( pa > 0_ni )
                read(val(2_ni:pa-1_ni), *) dshape(n)
                val(1_ni:pa) = ' '
                pa = index(val, ',')
                n = n + 1_nb
             end do
             read(val(2_ni:),*) dshape(n)
          case default
             write(*,*) 'Warning: Ignoring unknown key', key, ' with value ', val
       end select
    end do
    !
  end subroutine
  !
  ! Allocate a 4-dimensonal array
  subroutine allocate_f8_4d_by_shape(dat, dshape)
    real(kind=nr), intent(out), allocatable :: dat(:,:,:,:)
    integer(kind=ni), intent(in) :: dshape(:) 
    ! -----------------------------------------------------------------
    !
    if ( size(dshape) /= 4_ni ) then
       write(*,*) 'Error allocating: given shape is not 4D: ', dshape
       return
    end if
    !
    ! Fortran at it's best.
    allocate(dat(dshape(1_ni),dshape(2_ni),dshape(3_ni),dshape(4_ni)))
    !
  end subroutine
  !
  ! Allocate a 2-dimensonal array
  subroutine allocate_f8_2d_by_shape(dat, dshape)
    real(kind=nr), intent(out), allocatable :: dat(:,:)
    integer(kind=ni), intent(in) :: dshape(:)
    ! -----------------------------------------------------------------
    !
    if ( size(dshape) /= 2_ni ) then
       write(*,*) 'Error allocating: given shape is not 2D: ', dshape
       return
    end if
    !
    ! Fortran at it's best.
    allocate(dat(dshape(1_ni),dshape(2_ni)))
    !
  end subroutine
  !
  ! Read a 4-dimensional 64-bit real array
  subroutine read_f8_4d(runit, dat)
    real(kind=nr), intent(out) :: dat(:,:,:,:)
    integer(kind=ni), intent(in) :: runit
    ! -----------------------------------------------------------------
    !
    read(runit, pos=startpos(runit-runit0+1_ni)) dat
    !
  end subroutine
  !
  ! Close a file; Currently the runit will nevertheless not be reused.
  subroutine close_npy(runit)
    integer(kind=ni), intent(in) :: runit
    ! -----------------------------------------------------------------
    !
    close(runit)
  end subroutine
  !
end module
