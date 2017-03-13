module MOM_transform_test

!***********************************************************************
!*                   GNU General Public License                        *
!* This file is a part of MOM.                                         *
!*                                                                     *
!* MOM is free software; you can redistribute it and/or modify it and  *
!* are expected to follow the terms of the GNU General Public License  *
!* as published by the Free Software Foundation; either version 2 of   *
!* the License, or (at your option) any later version.                 *
!*                                                                     *
!* MOM is distributed in the hope that it will be useful, but WITHOUT  *
!* ANY WARRANTY; without even the implied warranty of MERCHANTABILITY  *
!* or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public    *
!* License for more details.                                           *
!*                                                                     *
!* For the full text of the GNU General Public License,                *
!* write to: Free Software Foundation, Inc.,                           *
!*           675 Mass Ave, Cambridge, MA 02139, USA.                   *
!* or see:   http://www.gnu.org/licenses/gpl.html                      *
!***********************************************************************

use MOM_coms, only : PE_here, root_PE
use MOM_error_handler, only : MOM_error, FATAL
use MOM_file_parser, only : log_version, get_param, param_file_type
use MOM_error_handler,  only : callTree_enter, callTree_leave

use mpp_mod, only : mpp_gather, mpp_max
use ensemble_manager_mod, only : get_ensemble_size, get_ensemble_id, get_ensemble_pelist

implicit none ; private

public :: MOM_transform_test_init, transform_test_started
public :: transform, transform_and_swap
public :: transform_allocatable, transform_allocatable_and_swap
public :: do_transform_test, do_transform_on_this_pe
public :: transform_compare, undo_transform

interface transform
  module procedure transform_2d, transform_3d
end interface

interface transform_allocatable
  module procedure transform_allocatable_2d, transform_allocatable_3d, &
                   transform_allocatable_4d
end interface

interface transform_allocatable_and_swap
  module procedure transform_allocatable_and_swap_1d, &
                   transform_allocatable_and_swap_2d, &
                   transform_allocatable_and_swap_3d
end interface

interface undo_transform
  module procedure undo_transform_2d, undo_transform_3d
end interface

interface transform_compare
  module procedure transform_compare_1d, transform_compare_2d, transform_compare_3d
end interface

interface transform_and_swap
  module procedure transform_and_swap_2d, transform_and_swap_3d
end interface

!> Whether or not we're in a transform test run
logical :: transform_test = .false.
!> Whether the transform being done on this PE?
logical :: transform_on_this_pe = .false.
!> Whether the test has started. No comparisons are done before this
! flag is set.
logical :: test_started = .false.

contains

! =====================================================================

!> MOM_transform_test_init initializes the module.
subroutine MOM_transform_test_init(param_file)
  type(param_file_type),   intent(in)    :: param_file
! This include declares and sets the variable "version".
#include "version_variable.h"
  character(len=40)  :: mod = "MOM_transform_test" ! This module's name.

  integer, dimension(6) :: ensemble_size
  integer, dimension(:, :), allocatable :: ensemble_pelist

  if (test_started) then
    return
  endif

  call log_version(param_file, mod, version)

  call get_param(param_file, mod, "TRANSFORM_TEST", &
                 transform_test, &
                 "Whether or not to run a transformation \n"//&
                 "test. This involves transposing or rotating to all \n"//&
                 "model inputs. This is a testing feature that can be \n"//&
                 "used to help find horizontal indexing errors. \n", default=.false.)

  ! Check that we're running as an ensemble. And that each ensemble member uses 1 PE.
  if (transform_test) then
    ensemble_size = get_ensemble_size()
    if (ensemble_size(1) == 1) then
      call MOM_error(FATAL, &
                     "TRANSFORM_TEST: Model not within ensemble")
    endif
    if (ensemble_size(2) /= 1) then
      call MOM_error(FATAL, &
                     "TRANSFORM_TEST: must have 1 PE per ensemble")
    endif

    allocate(ensemble_pelist(ensemble_size(1), ensemble_size(2)))
    call get_ensemble_pelist(ensemble_pelist)

    ! For this test the root PE will be the transformed run and the other
    ! will be the vanilla run.
    if (PE_here() == ensemble_pelist(1, 1)) then
      transform_on_this_pe = .true.
    endif

    deallocate(ensemble_pelist)

  endif

  test_started = .true.

end subroutine MOM_transform_test_init

function transform_test_started()
    logical :: transform_test_started

    transform_test_started = test_started

end function transform_test_started

function do_transform_on_this_pe()
    logical :: do_transform_on_this_pe

    do_transform_on_this_pe = transform_on_this_pe

end function do_transform_on_this_pe

function do_transform_test()
    logical :: do_transform_test

    do_transform_test = transform_test

end function do_transform_test

!< 
subroutine transform_2d(arrayIn, arrayOut)
  real, dimension(:,:), intent(in) :: arrayIn !< Input array
  real, dimension(:,:), intent(out) :: arrayOut !< Transformed input

  if (.not. transform_on_this_pe) then
    call MOM_error(FATAL, 'transform_2d: should not be called on this PE.')
  endif

  if (.true.) then
      call transpose_2d(arrayIn, arrayOut)
  else
      ! Try a 90 degree rotation
      call rot90_2d(arrayIn, arrayOut, 1)
  endif

end subroutine transform_2d

subroutine transform_3d(arrayIn, arrayOut)
  real, dimension(:,:,:), intent(in) :: arrayIn !< Input array
  real, dimension(:,:,:), intent(out) :: arrayOut !< Transformed input

  if (.not. transform_on_this_pe) then
    call MOM_error(FATAL, 'transform_3d: should not be called on this PE.')
  endif

  ! Try a 90 degree rotation
  if (.true.) then
      call transpose_3d(arrayIn, arrayOut)
  else
      call rot90_3d(arrayIn, arrayOut, 1)
  endif

end subroutine transform_3d

!< Transform an allocatable array. After this call input may have
! a different shape.
subroutine transform_allocatable_2d(array)
  real, dimension(:,:), allocatable, intent(inout) :: array

  real, allocatable, dimension(:,:) :: tmp
  integer :: isz, jsz

  if (.not. transform_on_this_pe) then
    call MOM_error(FATAL, 'transform_2d: should not be called on this PE.')
  endif

  isz = size(array, 1)
  jsz = size(array, 2)

  allocate(tmp(isz, jsz))
  tmp(:, :) = array(:, :)
  deallocate(array)
  allocate(array(jsz, isz))

  if (.true.) then
      call transpose_2d(tmp, array)
  else
      ! Try a 90 degree rotation
      call rot90_2d(tmp, array, 1)
  endif

  deallocate(tmp)

end subroutine transform_allocatable_2d

!< Transform an allocatable array. After this call input may have
! a different shape.
subroutine transform_allocatable_3d(array)
  real, dimension(:, :, :), allocatable, intent(inout) :: array

  real, allocatable, dimension(:, :, :) :: tmp
  integer :: isz, jsz, ksz

  if (.not. transform_on_this_pe) then
    call MOM_error(FATAL, 'transform_3d: should not be called on this PE.')
  endif

  isz = size(array, 1)
  jsz = size(array, 2)
  ksz = size(array, 3)

  allocate(tmp(isz, jsz, ksz))
  tmp(:, :, :) = array(:, :, :)
  deallocate(array)
  allocate(array(jsz, isz, ksz))

  if (.true.) then
      call transpose_3d(tmp, array)
  else
      ! Try a 90 degree rotation
      call rot90_3d(tmp, array, 1)
  endif

  deallocate(tmp)

end subroutine transform_allocatable_3d

!< Transform an allocatable array. After this call input may have
! a different shape.
subroutine transform_allocatable_4d(array)
  real, dimension(:, :, :, :), allocatable, intent(inout) :: array

  real, allocatable, dimension(:, :, :, :) :: tmp
  integer :: isz, jsz, ksz, lsz

  if (.not. transform_on_this_pe) then
    call MOM_error(FATAL, 'transform_4d: should not be called on this PE.')
  endif

  isz = size(array, 1)
  jsz = size(array, 2)
  ksz = size(array, 3)
  lsz = size(array, 4)

  allocate(tmp(isz, jsz, ksz, lsz))
  tmp(:, :, :, :) = array(:, :, :, :)
  deallocate(array)
  allocate(array(jsz, isz, ksz, lsz))

  if (.true.) then
      call transpose_4d(tmp, array)
  else
      ! Try a 90 degree rotation
      call MOM_error(FATAL, 'transform_4d: NotImplemented.')
  endif

  deallocate(tmp)

end subroutine transform_allocatable_4d

!< Transform an allocatable arrays and swap contents.
! After this call input may have a different shape.
subroutine transform_allocatable_and_swap_3d(arrayA, arrayB)
  real, dimension(:,:,:), allocatable, intent(inout) :: arrayA
  real, dimension(:,:,:), allocatable, intent(inout) :: arrayB

  real, allocatable, dimension(:,:,:) :: tmp
  integer :: isz, jsz, ksz

  if (.not. transform_on_this_pe) then
    call MOM_error(FATAL, 'transform_allocatable_and_swap_3d: should not be called on this PE.')
  endif

  if (size(arrayA, 1) /= size(arrayB, 1) .or. \
      size(arrayA, 2) /= size(arrayB, 2)) then
    call MOM_error(FATAL, 'transform_allocatable_and_swap_3d: array shapes  not compatible')
  endif

  isz = size(arrayA, 1)
  jsz = size(arrayA, 2)
  ksz = size(arrayA, 3)

  allocate(tmp(isz, jsz, ksz))
  tmp(:, :, :) = arrayA(:, :, :)
  deallocate(arrayA)
  allocate(arrayA(jsz, isz, ksz))

  if (.true.) then
      call transpose_3d(arrayB, arrayA)
  else
      ! Try a 90 degree rotation
      call rot90_3d(tmp, arrayA, 1)
  endif

  deallocate(arrayB)
  allocate(arrayB(jsz, isz, ksz))

  if (.true.) then
      call transpose_3d(tmp, arrayB)
  else
      ! Try a 90 degree rotation
      call rot90_3d(tmp, arrayB, 1)
  endif

  deallocate(tmp)

end subroutine transform_allocatable_and_swap_3d

!< Transform an allocatable arrays and swap contents.
! After this call input may have a different shape.
subroutine transform_allocatable_and_swap_2d(arrayA, arrayB)
  real, dimension(:,:), allocatable, intent(inout) :: arrayA
  real, dimension(:,:), allocatable, intent(inout) :: arrayB

  real, allocatable, dimension(:,:) :: tmp
  integer :: isz, jsz

  if (.not. transform_on_this_pe) then
    call MOM_error(FATAL, 'transform_allocatable_and_swap_2d: should not be called on this PE.')
  endif

  if (size(arrayA, 1) /= size(arrayB, 1) .or. \
      size(arrayA, 2) /= size(arrayB, 2)) then
    call MOM_error(FATAL, 'transform_allocatable_and_swap_2d: array shapes  not compatible')
  endif

  isz = size(arrayA, 1)
  jsz = size(arrayA, 2)

  allocate(tmp(isz, jsz))
  tmp(:, :) = arrayA(:, :)
  deallocate(arrayA)
  allocate(arrayA(jsz, isz))

  if (.true.) then
      call transpose_2d(arrayB, arrayA)
  else
      ! Try a 90 degree rotation
      call rot90_2d(tmp, arrayA, 1)
  endif

  deallocate(arrayB)
  allocate(arrayB(jsz, isz))

  if (.true.) then
      call transpose_2d(tmp, arrayB)
  else
      ! Try a 90 degree rotation
      call rot90_2d(tmp, arrayB, 1)
  endif

  deallocate(tmp)

end subroutine transform_allocatable_and_swap_2d

!< There is no 1d transformation, so this just swaps.
subroutine transform_allocatable_and_swap_1d(arrayA, arrayB)
  real, allocatable, dimension(:), intent(inout) :: arrayA, arrayB

  real, allocatable, dimension(:) :: tmp

  if (.not. transform_on_this_pe) then
    call MOM_error(FATAL, 'transform_allocatable_and_swap_1d: should not be called on this PE.')
  endif

  allocate(tmp(size(arrayA)))

  tmp(:) = arrayA(:)

  deallocate(arrayA)
  allocate(arrayA(size(arrayB)))
  arrayA(:) = arrayB

  deallocate(arrayB)
  allocate(arrayB(size(tmp)))
  arrayB(:) = tmp(:)

  deallocate(tmp)

end subroutine transform_allocatable_and_swap_1d


subroutine undo_transform_2d(original, undone)
  real, dimension(:,:), intent(in) :: original  !< The transformed array
  real, dimension(:,:), intent(out) :: undone !< The un-transformed array

  if (.not. transform_on_this_pe) then
    call MOM_error(FATAL, 'undo_transform_2d: should not be called on this PE.')
  endif

  if (.true.) then
      call transpose_2d(original, undone)
  else
      ! Try a 90 degree rotation
      call rot90_2d(original, undone, 3)
  endif

end subroutine undo_transform_2d


subroutine undo_transform_3d(original, reversed)
  real, dimension(:,:,:), intent(in) :: original  !< The transformed array
  real, dimension(:,:,:), intent(out) :: reversed !< The un-transformed array

  if (.not. transform_on_this_pe) then
    call MOM_error(FATAL, 'undo_transform_3d: should not be called on this PE.')
  endif

  if (.true.) then
      call transpose_3d(original, reversed)
  else
      ! Try a 90 degree rotation
      call rot90_3d(original, reversed, 3)
  endif

end subroutine undo_transform_3d

subroutine transpose_2d(arrayIn, arrayOut)
  real, dimension(:,:), intent(in) :: arrayIn !< Array to be transposed
  real, dimension(:,:), intent(out) :: arrayOut !< Transposed array

  if (size(arrayIn, 1) /= size(arrayOut, 2) &
      .or. size(arrayIn, 2) /= size(arrayOut, 1)) then
    call MOM_error(FATAL, 'transform_2d: array shapes incompatible.')
  endif

  arrayOut(:, :) = transpose(arrayIn(:, :))

end subroutine transpose_2d

subroutine transpose_3d(arrayIn, arrayOut)
  real, dimension(:,:,:), intent(in) :: arrayIn !< The array to be transposed
  real, dimension(:,:,:), intent(inout) :: arrayOut !< Transposed array

  integer :: k

  do k=lbound(arrayIn, 3), ubound(arrayIn, 3)
     call transpose_2d(arrayIn(:, :, k), arrayOut(:, :, k))
  enddo

end subroutine transpose_3d

subroutine transpose_4d(arrayIn, arrayOut)
  real, dimension(:,:,:,:), intent(in) :: arrayIn !< The array to be transposed
  real, dimension(:,:,:,:), intent(inout) :: arrayOut !< Transposed array

  integer :: l

  do l=lbound(arrayIn, 4), ubound(arrayIn, 4)
     call transpose_3d(arrayIn(:, :, : , l), arrayOut(:, :, :, l))
  enddo

end subroutine transpose_4d

subroutine transform_and_swap_2d(arrayA, arrayB)
  real, intent(inout), dimension(:,:) :: arrayA, arrayB

  real, allocatable, dimension(:,:) :: tmp

  if (.not. transform_on_this_pe) then
    call MOM_error(FATAL, 'transform_and_swap_2d: should not be called on this PE.')
  endif

  if (size(arrayA, 1) /= size(arrayB, 1) .or. size(arrayA, 2) /= size(arrayB, 2)) then
    call MOM_error(FATAL, 'transform_and_swap_2d: arrays shapes incompatible.')
  endif

  allocate(tmp(size(arrayA, 1), size(arrayA, 2)))

  tmp(:, :) = arrayA(:, :)

  if (.true.) then
    call transpose_2d(arrayB, arrayA)
    call transpose_2d(tmp, arrayB)
  else
    call rot90_2d(arrayB, arrayA, 1)
    call rot90_2d(tmp, arrayB, 1)
  endif

  deallocate(tmp)

end subroutine transform_and_swap_2d

subroutine transform_and_swap_3d(arrayA, arrayB)
  real, intent(inout), dimension(:,:,:) :: arrayA, arrayB

  real, allocatable, dimension(:,:,:) :: tmp

  if (.not. transform_on_this_pe) then
    call MOM_error(FATAL, 'transform_and_swap_3d: should not be called on this PE.')
  endif

  if (size(arrayA, 1) /= size(arrayB, 2) .or. &
          size(arrayA, 2) /= size(arrayB, 1) .or. &
          size(arrayA, 3) /= size(arrayB, 3)) then
    call MOM_error(FATAL, 'trans_and_swap_3d: array shapes incompatible.')
  endif

  allocate(tmp(size(arrayA, 1), size(arrayA, 2), size(arrayA, 3)))

  tmp(:, :, :) = arrayA(:, :, :)

  if (.true.) then
    call transpose_3d(arrayB, arrayA)
    call transpose_3d(tmp, arrayB)
  else
    call rot90_3d(arrayB, arrayA, 1)
    call rot90_3d(tmp, arrayB, 1)
  endif

  deallocate(tmp)

end subroutine transform_and_swap_3d

subroutine transform_compare_2d(arrayA, arrayB, ret)
  real, dimension(:, :), intent(in) :: arrayA, arrayB
  integer, intent(out) :: ret

  real, allocatable, dimension(:,:) :: tmp

  ret = 1
  if (.not. test_started) then
    ret = 0
    return
  endif

  if (transform_on_this_pe) then
    allocate(tmp(size(arrayA, 2), size(arrayA, 1)))
    call undo_transform_2d(arrayA, tmp)
    call ensemble_compare_1d(reshape(tmp, (/ size(tmp) /)), ret)

    if (ret /= 0) then
      call write_to_netcdf_2d(tmp, 'transform_test_debug_A.nc')
      deallocate(tmp)
      return
    endif

    deallocate(tmp)

    allocate(tmp(size(arrayB, 2), size(arrayB, 1)))
    call undo_transform_2d(arrayB, tmp)
    call ensemble_compare_1d(reshape(tmp, (/ size(tmp) /)), ret)

    if (ret /= 0) then
      call write_to_netcdf_2d(tmp, 'transform_test_debug_B.nc')
    endif
    deallocate(tmp)
  else
    call ensemble_compare_1d(reshape(arrayB, (/ size(arrayB) /)), ret)

    if (ret /= 0) then
      call write_to_netcdf_2d(arrayB, 'transform_test_debug_B.nc')
      return
    endif

    call ensemble_compare_1d(reshape(arrayA, (/ size(arrayA) /)), ret)

    if (ret /= 0) then
      call write_to_netcdf_2d(arrayB, 'transform_test_debug_A.nc')
    endif
  endif

end subroutine transform_compare_2d

subroutine transform_compare_3d(arrayA, arrayB, ret)
  real, dimension(:, :, :), intent(in) :: arrayA, arrayB
  integer, intent(out) :: ret

  real, allocatable, dimension(:, :, :) :: tmp

  if (.not. test_started) then
    ret = 0
    return
  endif

  if (transform_on_this_pe) then
    allocate(tmp(size(arrayA, 2), size(arrayA, 1), size(arrayA, 3)))
    call undo_transform_3d(arrayA, tmp)
    call ensemble_compare_1d(reshape(tmp, (/ size(tmp) /)), ret)

    if (ret /= 0) then
      call write_to_netcdf_3d(tmp, 'transform_test_debug_A.nc')
    endif
    deallocate(tmp)
  else
    call ensemble_compare_1d(reshape(arrayB, (/ size(arrayB) /)), ret)

    if (ret /= 0) then
      call write_to_netcdf_3d(arrayB, 'transform_test_debug_B.nc')
    endif
  endif

end subroutine transform_compare_3d

subroutine transform_compare_1d(arrayA, arrayB, ret)
  real, intent(in), dimension(:) :: arrayA, arrayB
  integer, intent(out) :: ret

  if (.not. test_started) then
    ret = 0
    return
  endif

  if (transform_on_this_pe) then
    call ensemble_compare_1d(arrayA, ret)
  else
    call ensemble_compare_1d(arrayB, ret)
  endif

end subroutine transform_compare_1d

subroutine ensemble_compare_1d(sbuf, ret)
  real, intent(in), dimension(:) :: sbuf
  integer, intent(out) :: ret

  integer, dimension(:, :), allocatable :: ensemble_pelist
  integer, dimension(6) :: ensemble_size
  real, dimension(:), allocatable :: rbuf
  integer :: sbuf_size, e, i
  real :: a, b

  ret = 0

  ! If we are running in an ensemble then communicate with other member.
  ensemble_size = get_ensemble_size()
  if (ensemble_size(1) == 1) then
    return
  endif

  ! ensemble_size(1) is the number of ensemble members
  ! ensemble_size(2) is the number of pes per ensemble
  sbuf_size = size(sbuf)
  allocate(rbuf(ensemble_size(1) * sbuf_size))

  allocate(ensemble_pelist(ensemble_size(1), ensemble_size(2)))
  call get_ensemble_pelist(ensemble_pelist)

  ! Gather to root of every ensemble member.
  call mpp_gather(sbuf, rbuf, ensemble_pelist(:, 1))

  ! Check that sbuf is the same on all pes.
  if (transform_on_this_pe) then
    do e=0,ensemble_size(1)-1
      do i=1,sbuf_size
        a = sbuf(i)
        b = rbuf(e*sbuf_size + i)
        if ((a /= b) .and. (abs(a - b) /= 0.0)) then
          print*, a - b
          ret = i
          exit
        endif
      enddo
    enddo
  endif

  deallocate(rbuf)

  call mpp_max(ret, ensemble_pelist(:, 1))

end subroutine ensemble_compare_1d

subroutine write_to_netcdf_3d(array, file_name)
  use netcdf
  implicit none

  real, intent(in), dimension(:, :,:) :: array
  character(len=*), intent(in) :: file_name

  integer :: file_id, xdim_id, ydim_id, zdim_id
  integer :: array_id
  integer, dimension(3) :: arrdims
  character(len=*), parameter :: arrunit = 'ergs'

  integer :: i, j, k
  integer :: ierr

  i = size(array,1)
  j = size(array,2)
  k = size(array,3)

  ! create the file
  ierr = nf90_create(path=trim(file_name), cmode=NF90_CLOBBER, ncid=file_id)

  ! define the dimensions
  ierr = nf90_def_dim(file_id, 'X', i, xdim_id)
  ierr = nf90_def_dim(file_id, 'Y', j, ydim_id)
  ierr = nf90_def_dim(file_id, 'Z', k, zdim_id)

  ! now that the dimensions are defined, we can define variables on them,...
  arrdims = (/ xdim_id, ydim_id, zdim_id /)
  ierr = nf90_def_var(file_id, 'Array',  NF90_DOUBLE, arrdims, array_id)

  ! ...and assign units to them as an attribute
  ierr = nf90_put_att(file_id, array_id, "units", arrunit)

  ! done defining
  ierr = nf90_enddef(file_id)

  ! Write out the values
  ierr = nf90_put_var(file_id, array_id, array)

  ! close; done
  ierr = nf90_close(file_id)
end subroutine write_to_netcdf_3d

subroutine write_to_netcdf_2d(array, file_name)
  use netcdf
  implicit none

  real, intent(in), dimension(:,:) :: array
  character(len=*), intent(in) :: file_name

  integer :: file_id, xdim_id, ydim_id
  integer :: array_id
  integer, dimension(2) :: arrdims
  character(len=*), parameter :: arrunit = 'ergs'

  integer :: i, j
  integer :: ierr

  i = size(array,1)
  j = size(array,2)

  ! create the file
  ierr = nf90_create(path=trim(file_name), cmode=NF90_CLOBBER, ncid=file_id)

  ! define the dimensions
  ierr = nf90_def_dim(file_id, 'X', i, xdim_id)
  ierr = nf90_def_dim(file_id, 'Y', j, ydim_id)

  ! now that the dimensions are defined, we can define variables on them,...
  arrdims = (/ xdim_id, ydim_id /)
  ierr = nf90_def_var(file_id, 'Array',  NF90_DOUBLE, arrdims, array_id)

  ! ...and assign units to them as an attribute
  ierr = nf90_put_att(file_id, array_id, "units", arrunit)

  ! done defining
  ierr = nf90_enddef(file_id)

  ! Write out the values
  ierr = nf90_put_var(file_id, array_id, array)

  ! close; done
  ierr = nf90_close(file_id)
end subroutine write_to_netcdf_2d

subroutine rot90_2d(arrayIn, arrayOut, nrot90)
  real, dimension(:,:), intent(in) :: arrayIn !< Array to be rotated
  real, dimension(:,:), intent(out) :: arrayOut !< Rotated array
  integer, intent(in) :: nrot90 !< Number of 90 degree rotations to perform

  if (.not. nrot90 < 4) then
    call MOM_error(FATAL, 'rot90_2d: nrot should be < 4')
  endif

  if (modulo(nrot90, 2) == 2) then
    if (size(arrayIn, 1) /= size(arrayOut, 1) .or. size(arrayIn, 2) /= size(arrayOut, 2)) then
      call MOM_error(FATAL, 'rot90_2d: 180 deg rotation bad array shapes.')
    endif
  else
    if (size(arrayIn, 1) /= size(arrayOut, 2) .or. size(arrayIn, 2) /= size(arrayOut, 1)) then
      call MOM_error(FATAL, 'rot90_2d: 90 deg rotation bad array shapes.')
    endif
  endif

  if (nrot90 == 1) then
    ! transpose, reverse rows
    arrayOut(:, :) = transpose(arrayIn(:, :))
    arrayOut(:, :) = arrayOut(:, ubound(arrayOut, 2):lbound(arrayOut, 2):-1)
  elseif (nrot90 == 2) then
    ! reverse both rows and cols
    arrayOut(:, :) = arrayIn(ubound(arrayIn, 1):lbound(arrayIn, 1):-1, &
                             ubound(arrayIn, 2):lbound(arrayIn, 2):-1)
  elseif (nrot90 == 3) then
    ! transpose, reverse cols
    arrayOut(:,:) = transpose(arrayIn(:, :))
    arrayOut(:, :) = arrayOut(ubound(arrayOut, 1):lbound(arrayOut, 1):-1, :)
  endif

end subroutine rot90_2d

subroutine rot90_3d(arrayIn, arrayOut, nrot90)
  real, dimension(:,:,:), intent(in) :: arrayIn !< The array to be rotated
  real, dimension(:,:,:), intent(inout) :: arrayOut !< Rotated array
  integer, intent(in) :: nrot90 !< Number of 90 degree rotations to perform

  integer :: k

  if (.not. nrot90 < 4) then
    call MOM_error(FATAL, 'rot90_2d: nrot should be < 4')
  endif

  do k=lbound(arrayIn, 3), ubound(arrayIn, 3)
     call rot90_2d(arrayIn(:, :, k), arrayOut(:, :, k), nrot90)
  enddo

end subroutine rot90_3d

end module MOM_transform_test
