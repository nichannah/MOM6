module MOM_checksums

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

use MOM_coms, only : PE_here, root_PE, num_PEs, sum_across_PEs
use MOM_coms, only : min_across_PEs, max_across_PEs
use MOM_coms, only : reproducing_sum
use MOM_error_handler, only : MOM_error, FATAL, is_root_pe
use MOM_file_parser, only : log_version, get_param, param_file_type
use MOM_hor_index, only : hor_index_type

use mpp_mod, only : mpp_gather
use ensemble_manager_mod, only : get_ensemble_size, get_ensemble_id, get_ensemble_pelist

implicit none ; private

public :: hchksum, Bchksum, uchksum, vchksum, qchksum, chksum, is_NaN
public :: hchksum_pair, uvchksum_pair, Bchksum_pair
public :: transform_input, transform_and_swap_input, do_transform_input
public :: MOM_checksums_init

interface hchksum_pair
  module procedure chksum_pair_h_2d, chksum_pair_h_3d
end interface

interface uvchksum_pair
  module procedure chksum_pair_uv_2d, chksum_pair_uv_3d
end interface

interface Bchksum_pair
  module procedure chksum_pair_B_2d, chksum_pair_B_3d
end interface

interface hchksum
  module procedure chksum_h_2d, chksum_h_3d
end interface

interface Bchksum
  module procedure chksum_B_2d, chksum_B_3d
end interface

interface qchksum
  module procedure chksum_B_2d, chksum_B_3d
end interface

interface uchksum
  module procedure chksum_u_2d, chksum_u_3d
end interface

interface vchksum
  module procedure chksum_v_2d, chksum_v_3d
end interface

interface chksum
  module procedure chksum1d, chksum2d, chksum3d
end interface

interface chk_sum_msg
  module procedure chk_sum_msg1, chk_sum_msg2, chk_sum_msg3, chk_sum_msg5
end interface

interface is_NaN
  module procedure is_NaN_0d, is_NaN_1d, is_NaN_2d, is_NaN_3d
end interface

interface transform_input
  module procedure transform_input_2d, transform_input_3d
end interface

interface transform_and_swap_input
  module procedure transform_and_swap_input_2d, transform_and_swap_input_3d
end interface

integer, parameter :: default_shift=0
logical :: calculateStatistics=.true. ! If true, report min, max and mean.
logical :: writeChksums=.true. ! If true, report the bitcount checksum
logical :: checkForNaNs=.true. ! If true, checks array for NaNs and cause
                               ! FATAL error is any are found
logical :: transform_input_test = .false.
logical :: transform_input_on_this_pe = .false.

contains

! =====================================================================

subroutine chksum_pair_h_2d(arrayA, mesgA, arrayB, mesgB, HI, haloshift)
  type(hor_index_type),             intent(in) :: HI     !< A horizontal index type
  real, dimension(HI%isd:,HI%jsd:), intent(in) :: arrayA, arrayB !< The arrays to be checksummed
  character(len=*),                 intent(in) :: mesgA, mesgB !< Identifying messages
  integer,                optional, intent(in) :: haloshift !< The width of halos to check (default 0)

  if (transform_input_test) then
    if (transform_input_on_this_pe) then
      call compare_within_ensemble(reshape(transpose(arrayA), (/ size(arrayA) /)))
    else
      call compare_within_ensemble(reshape(arrayB, (/ size(arrayB) /)))
    endif
  endif

  if (present(haloshift)) then
    call chksum_h_2d(arrayA, mesgA, HI, haloshift, compare=.false.)
    call chksum_h_2d(arrayB, mesgB, HI, haloshift, compare=.false.)
  else
    call chksum_h_2d(arrayA, mesgA, HI, compare=.false.)
    call chksum_h_2d(arrayB, mesgB, HI, compare=.false.)
  endif

end subroutine chksum_pair_h_2d

subroutine chksum_pair_h_3d(arrayA, mesgA, arrayB, mesgB, HI, haloshift)
  type(hor_index_type),                intent(in) :: HI     !< A horizontal index type
  real, dimension(HI%isd:,HI%jsd:, :), intent(in) :: arrayA, arrayB !< The arrays to be checksummed
  character(len=*),                    intent(in) :: mesgA, mesgB !< Identifying messages
  integer,                   optional, intent(in) :: haloshift !< The width of halos to check (default 0)

end subroutine chksum_pair_h_3d

!> chksum_h_2d performs checksums on a 2d array staggered at tracer points.
subroutine chksum_h_2d(array, mesg, HI, haloshift, fname, compare)
  type(hor_index_type),           intent(in) :: HI     !< A horizontal index type
  real, dimension(HI%isd:,HI%jsd:), intent(in) :: array !< The array to be checksummed
  character(len=*),                intent(in) :: mesg  !< An identifying message
  integer,               optional, intent(in) :: haloshift !< The width of halos to check (default 0)
  character(len=*),      optional, intent(in) :: fname  !< Name of file to dump
  logical,               optional, intent(in) :: compare  !< Compare if in the transform input test

  integer :: bc0,bcSW,bcSE,bcNW,bcNE,hshift
  logical :: do_compare

  do_compare = .true.
  if (present(compare)) then
    do_compare = compare
  endif

  if (present(fname) .and. transform_input_test) then
      if (transform_input_on_this_pe) then
        call write_to_netcdf_2d(array, fname//'_h_trans.nc')
      else
        call write_to_netcdf_2d(array, fname//'_h.nc')
      endif
  endif

  if (do_compare .and. transform_input_test) then
    if (transform_input_on_this_pe) then
      call compare_within_ensemble(reshape(transpose(array), (/ size(array) /)))
    else
      call compare_within_ensemble(reshape(array, (/ size(array) /)))
    endif
  endif

  if (checkForNaNs) then
    if (is_NaN(array(HI%isc:HI%iec,HI%jsc:HI%jec))) &
      call chksum_error(FATAL, 'NaN detected: '//trim(mesg))
!   if (is_NaN(array)) &
!     call chksum_error(FATAL, 'NaN detected in halo: '//trim(mesg))
  endif

  if (calculateStatistics) call subStats(HI, array, mesg)

  if (.not.writeChksums) return

  hshift=default_shift
  if (present(haloshift)) hshift=haloshift
  if (hshift<0) hshift=HI%ied-HI%iec

  if ( HI%isc-hshift<HI%isd .or. HI%iec+hshift>HI%ied .or. &
       HI%jsc-hshift<HI%jsd .or. HI%jec+hshift>HI%jed ) then
    write(0,*) 'chksum_h_2d: haloshift =',hshift
    write(0,*) 'chksum_h_2d: isd,isc,iec,ied=',HI%isd,HI%isc,HI%iec,HI%ied
    write(0,*) 'chksum_h_2d: jsd,jsc,jec,jed=',HI%jsd,HI%jsc,HI%jec,HI%jed
    call chksum_error(FATAL,'Error in chksum_h_2d '//trim(mesg))
  endif

  bc0=subchk(array, HI, 0, 0)

  if (hshift==0) then
      if (is_root_pe()) call chk_sum_msg("h-point:",bc0,mesg)
      return
  endif

  bcSW=subchk(array, HI, -hshift, -hshift)
  bcSE=subchk(array, HI, hshift, -hshift)
  bcNW=subchk(array, HI, -hshift, hshift)
  bcNE=subchk(array, HI, hshift, hshift)

  if (is_root_pe()) call chk_sum_msg("h-point:",bc0,bcSW,bcSE,bcNW,bcNE,mesg)

  contains

  integer function subchk(array, HI, di, dj)
    type(hor_index_type), intent(in) :: HI
    real, dimension(HI%isd:,HI%jsd:), intent(in) :: array
    integer, intent(in) :: di, dj
    integer :: bitcount, i, j, bc
    subchk = 0
    do j=HI%jsc+dj,HI%jec+dj; do i=HI%isc+di,HI%iec+di
        bc = bitcount(abs(array(i,j)))
        subchk = subchk + bc
    enddo; enddo
    call sum_across_PEs(subchk)
    subchk=mod(subchk,1000000000)
  end function subchk

  subroutine subStats(HI, array, mesg)
    type(hor_index_type), intent(in) :: HI
    real, dimension(HI%isd:,HI%jsd:), intent(in) :: array
    character(len=*), intent(in) :: mesg
    integer :: i, j, n
    real :: aMean, aMin, aMax

    aMin = array(HI%isc,HI%jsc)
    aMax = array(HI%isc,HI%jsc)
    n = 0
    do j=HI%jsc,HI%jec ; do i=HI%isc,HI%iec
      aMin = min(aMin, array(i,j))
      aMax = max(aMax, array(i,j))
      n = n + 1
    enddo ; enddo

    aMean = reproducing_sum(array(HI%isc:HI%iec,HI%jsc:HI%jec))
    call sum_across_PEs(n)
    call min_across_PEs(aMin)
    call max_across_PEs(aMax)
    aMean = aMean / real(n)
    if (is_root_pe()) call chk_sum_msg("h-point:",aMean,aMin,aMax,mesg)
  end subroutine subStats

end subroutine chksum_h_2d

! =====================================================================

subroutine chksum_pair_B_2d(arrayA, mesgA, arrayB, mesgB, HI, haloshift)
  type(hor_index_type),             intent(in) :: HI     !< A horizontal index type
  real, dimension(HI%isd:,HI%jsd:), intent(in) :: arrayA, arrayB !< The arrays to be checksummed
  character(len=*),                 intent(in) :: mesgA, mesgB !< Identifying messages
  integer,                optional, intent(in) :: haloshift !< The width of halos to check (default 0)

  if (transform_input_test) then
    if (transform_input_on_this_pe) then
      call compare_within_ensemble(reshape(transpose(arrayA), (/ size(arrayA) /)))
    else
      call compare_within_ensemble(reshape(arrayB, (/ size(arrayB) /)))
    endif
  endif

  if (present(haloshift)) then
    call chksum_u_2d(arrayA, mesgA, HI, haloshift)
    call chksum_v_2d(arrayB, mesgB, HI, haloshift)
  else
    call chksum_u_2d(arrayA, mesgA, HI)
    call chksum_v_2d(arrayB, mesgB, HI)
  endif

end subroutine chksum_pair_B_2d

subroutine chksum_pair_B_3d(arrayA, mesgA, arrayB, mesgB, HI, haloshift)
  type(hor_index_type),                intent(in) :: HI     !< A horizontal index type
  real, dimension(HI%isd:,HI%jsd:, :), intent(in) :: arrayA, arrayB !< The arrays to be checksummed
  character(len=*),                    intent(in) :: mesgA, mesgB !< Identifying messages
  integer,                   optional, intent(in) :: haloshift !< The width of halos to check (default 0)

end subroutine chksum_pair_B_3d

!> chksum_B_2d performs checksums on a 2d array staggered at corner points.
subroutine chksum_B_2d(array, mesg, HI, haloshift, symmetric, fname)
  type(hor_index_type), intent(in) :: HI     !< A horizontal index type
  real, dimension(HI%IsdB:,HI%JsdB:), &
                        intent(in) :: array !< The array to be checksummed
  character(len=*),     intent(in) :: mesg  !< An identifying message
  integer,    optional, intent(in) :: haloshift !< The width of halos to check (default 0)
  logical,    optional, intent(in) :: symmetric !< If true, do the checksums on the
                                                !! full symmetric computational domain.
  character(len=*),     optional, intent(in) :: fname  !< Name of file to dump

  integer :: bc0,bcSW,bcSE,bcNW,bcNE,hshift
  logical :: sym

  if (present(fname) .and. transform_input_test) then
      if (transform_input_on_this_pe) then
        call write_to_netcdf_2d(array, fname//'_B_trans.nc')
      else
        call write_to_netcdf_2d(array, fname//'_B.nc')
      endif
  endif

  if (checkForNaNs) then
    if (is_NaN(array(HI%IscB:HI%IecB,HI%JscB:HI%JecB))) &
      call chksum_error(FATAL, 'NaN detected: '//trim(mesg))
!   if (is_NaN(array)) &
!     call chksum_error(FATAL, 'NaN detected in halo: '//trim(mesg))
  endif

  if (calculateStatistics) call subStats(HI, array, mesg)

  if (.not.writeChksums) return

  hshift=default_shift
  if (present(haloshift)) hshift=haloshift
  if (hshift<0) hshift=HI%ied-HI%iec

  if ( HI%iscB-hshift<HI%isdB .or. HI%iecB+hshift>HI%iedB .or. &
       HI%jscB-hshift<HI%jsdB .or. HI%jecB+hshift>HI%jedB ) then
    write(0,*) 'chksum_B_2d: haloshift =',hshift
    write(0,*) 'chksum_B_2d: isd,isc,iec,ied=',HI%isdB,HI%iscB,HI%iecB,HI%iedB
    write(0,*) 'chksum_B_2d: jsd,jsc,jec,jed=',HI%jsdB,HI%jscB,HI%jecB,HI%jedB
    call chksum_error(FATAL,'Error in chksum_B_2d '//trim(mesg))
  endif

  sym = .false. ; if (present(symmetric)) sym = symmetric

  bc0=subchk(array, HI, 0, 0)

  if ((hshift==0) .and. .not.sym) then
    if (is_root_pe()) call chk_sum_msg("B-point:",bc0,mesg)
    return
  endif

  if (sym) then
    bcSW=subchk(array, HI, -hshift-1, -hshift-1)
    bcSE=subchk(array, HI, hshift, -hshift-1)
    bcNW=subchk(array, HI, -hshift-1, hshift)
  else
    bcSW=subchk(array, HI, -hshift, -hshift)
    bcSE=subchk(array, HI, hshift, -hshift)
    bcNW=subchk(array, HI, -hshift, hshift)
  endif
  bcNE=subchk(array, HI, hshift, hshift)

  if (is_root_pe()) call chk_sum_msg("B-point:",bc0,bcSW,bcSE,bcNW,bcNE,mesg)

  contains

  integer function subchk(array, HI, di, dj)
    type(hor_index_type), intent(in) :: HI
    real, dimension(HI%IsdB:,HI%JsdB:), intent(in) :: array
    integer, intent(in) :: di, dj
    integer :: bitcount, i, j, bc
    subchk = 0
    do j=HI%jscB+dj,HI%jecB+dj; do i=HI%iscB+di,HI%iecB+di
        bc = bitcount(abs(array(i,j)))
        subchk = subchk + bc
    enddo; enddo
    call sum_across_PEs(subchk)
    subchk=mod(subchk,1000000000)
  end function subchk

  subroutine subStats(HI, array, mesg)
    type(hor_index_type), intent(in) :: HI
    real, dimension(HI%IsdB:,HI%JsdB:), intent(in) :: array
    character(len=*), intent(in) :: mesg
    integer :: i, j, n
    real :: aMean, aMin, aMax
    aMean = 0.
    aMin = array(HI%iscB,HI%jscB)
    aMax = array(HI%iscB,HI%jscB)
    n = 0
    do j=HI%jscB,HI%jecB ; do i=HI%iscB,HI%iecB
      aMin = min(aMin, array(i,j))
      aMax = max(aMax, array(i,j))
      n = n + 1
    enddo ; enddo
    aMean = reproducing_sum(array(HI%iscB:HI%iecB,HI%jscB:HI%jecB))
    call sum_across_PEs(n)
    call min_across_PEs(aMin)
    call max_across_PEs(aMax)
    aMean = aMean / real(n)
    if (is_root_pe()) call chk_sum_msg("B-point:",aMean,aMin,aMax,mesg)
  end subroutine subStats

end subroutine chksum_B_2d

! =====================================================================

subroutine chksum_pair_uv_2d(arrayU, mesgU, arrayV, mesgV, HI, haloshift)
  type(hor_index_type),             intent(in) :: HI     !< A horizontal index type
  real, dimension(HI%isd:,HI%jsd:), intent(in) :: arrayU, arrayV !< The arrays to be checksummed
  character(len=*),                 intent(in) :: mesgU, mesgV !< Identifying messages
  integer,                optional, intent(in) :: haloshift !< The width of halos to check (default 0)

  if (transform_input_test) then
    if (transform_input_on_this_pe) then
      call compare_within_ensemble(reshape(transpose(arrayU), (/ size(arrayV) /)))
    else
      call compare_within_ensemble(reshape(arrayV, (/ size(arrayV) /)))
    endif
  endif

  if (present(haloshift)) then
    call chksum_u_2d(arrayU, mesgU, HI, haloshift)
    call chksum_v_2d(arrayV, mesgV, HI, haloshift)
  else
    call chksum_u_2d(arrayU, mesgU, HI)
    call chksum_v_2d(arrayV, mesgV, HI)
  endif

end subroutine chksum_pair_uv_2d

subroutine chksum_pair_uv_3d(arrayU, mesgU, arrayV, mesgV, HI, haloshift)
  type(hor_index_type),                intent(in) :: HI     !< A horizontal index type
  real, dimension(HI%isd:,HI%jsd:, :), intent(in) :: arrayU, arrayV !< The arrays to be checksummed
  character(len=*),                    intent(in) :: mesgU, mesgV !< Identifying messages
  integer,                   optional, intent(in) :: haloshift !< The width of halos to check (default 0)

end subroutine chksum_pair_uv_3d

!> chksum_u_2d performs checksums on a 2d array staggered at C-grid u points.
subroutine chksum_u_2d(array, mesg, HI, haloshift, fname)
  type(hor_index_type),           intent(in) :: HI     !< A horizontal index type
  real, dimension(HI%IsdB:,HI%jsd:), intent(in) :: array !< The array to be checksummed
  character(len=*),                intent(in) :: mesg  !< An identifying message
  integer,               optional, intent(in) :: haloshift !< The width of halos to check (default 0)
  character(len=*),        optional, intent(in) :: fname  !< Name of file to dump

  integer :: bc0,bcSW,bcSE,bcNW,bcNE,hshift

  if (present(fname) .and. transform_input_test) then
      if (transform_input_on_this_pe) then
        call write_to_netcdf_2d(array, fname//'_u_trans.nc')
      else
        call write_to_netcdf_2d(array, fname//'_u.nc')
      endif
  endif

  if (checkForNaNs) then
    if (is_NaN(array(HI%IscB:HI%IecB,HI%jsc:HI%jec))) &
      call chksum_error(FATAL, 'NaN detected: '//trim(mesg))
!   if (is_NaN(array)) &
!     call chksum_error(FATAL, 'NaN detected in halo: '//trim(mesg))
  endif

  if (calculateStatistics) call subStats(HI, array, mesg)

  if (.not.writeChksums) return

  hshift=default_shift
  if (present(haloshift)) hshift=haloshift
  if (hshift<0) hshift=HI%iedB-HI%iecB

  if ( HI%iscB-hshift<HI%isdB .or. HI%iecB+hshift>HI%iedB .or. &
       HI%jsc-hshift<HI%jsd .or. HI%jec+hshift>HI%jed ) then
    write(0,*) 'chksum_u_2d: haloshift =',hshift
    write(0,*) 'chksum_u_2d: isd,isc,iec,ied=',HI%isdB,HI%iscB,HI%iecB,HI%iedB
    write(0,*) 'chksum_u_2d: jsd,jsc,jec,jed=',HI%jsd,HI%jsc,HI%jec,HI%jed
    call chksum_error(FATAL,'Error in chksum_u_2d '//trim(mesg))
  endif

  bc0=subchk(array, HI, 0, 0)

  if (hshift==0) then
      if (is_root_pe()) call chk_sum_msg("u-point:",bc0,mesg)
      return
  endif

  bcSW=subchk(array, HI, -hshift, -hshift)
  bcSE=subchk(array, HI, hshift, -hshift)
  bcNW=subchk(array, HI, -hshift, hshift)
  bcNE=subchk(array, HI, hshift, hshift)

  if (is_root_pe()) call chk_sum_msg("u-point:",bc0,bcSW,bcSE,bcNW,bcNE,mesg)

  contains

  integer function subchk(array, HI, di, dj)
    type(hor_index_type), intent(in) :: HI
    real, dimension(HI%IsdB:,HI%jsd:), intent(in) :: array
    integer, intent(in) :: di, dj
    integer :: bitcount, i, j, bc
    subchk = 0
    do j=HI%jsc+dj,HI%jec+dj; do i=HI%iscB+di,HI%iecB+di
        bc = bitcount(abs(array(i,j)))
        subchk = subchk + bc
    enddo; enddo
    call sum_across_PEs(subchk)
    subchk=mod(subchk,1000000000)
  end function subchk

  subroutine subStats(HI, array, mesg)
    type(hor_index_type), intent(in) :: HI
    real, dimension(HI%IsdB:,HI%jsd:), intent(in) :: array
    character(len=*), intent(in) :: mesg
    integer :: i, j, n
    real :: aMean, aMin, aMax

    aMin = array(HI%iscB,HI%jsc)
    aMax = array(HI%iscB,HI%jsc)
    n = 0
    do j=HI%jsc,HI%jec ; do i=HI%iscB,HI%iecB
      aMin = min(aMin, array(i,j))
      aMax = max(aMax, array(i,j))
      n = n + 1
    enddo ; enddo
    aMean = reproducing_sum(array(HI%iscB:HI%iecB,HI%jsc:HI%jec))
    call sum_across_PEs(n)
    call min_across_PEs(aMin)
    call max_across_PEs(aMax)
    aMean = aMean / real(n)
    if (is_root_pe()) call chk_sum_msg("u-point:",aMean,aMin,aMax,mesg)
  end subroutine subStats

end subroutine chksum_u_2d

! =====================================================================

!> chksum_v_2d performs checksums on a 2d array staggered at C-grid v points.
subroutine chksum_v_2d(array, mesg, HI, haloshift, fname)
  type(hor_index_type),           intent(in) :: HI     !< A horizontal index type
  real, dimension(HI%isd:,HI%JsdB:), intent(in) :: array !< The array to be checksummed
  character(len=*),                intent(in) :: mesg  !< An identifying message
  integer,               optional, intent(in) :: haloshift !< The width of halos to check (default 0)
  character(len=*),        optional, intent(in) :: fname  !< Name of file to dump

  integer :: bc0,bcSW,bcSE,bcNW,bcNE,hshift

  if (present(fname) .and. transform_input_test) then
      if (transform_input_on_this_pe) then
        call write_to_netcdf_2d(array, fname//'_v_trans.nc')
      else
        call write_to_netcdf_2d(array, fname//'_v.nc')
      endif
  endif

  if (checkForNaNs) then
    if (is_NaN(array(HI%isc:HI%iec,HI%JscB:HI%JecB))) &
      call chksum_error(FATAL, 'NaN detected: '//trim(mesg))
!   if (is_NaN(array)) &
!     call chksum_error(FATAL, 'NaN detected in halo: '//trim(mesg))
  endif

  if (calculateStatistics) call subStats(HI, array, mesg)

  if (.not.writeChksums) return

  hshift=default_shift
  if (present(haloshift)) hshift=haloshift
  if (hshift<0) hshift=HI%ied-HI%iec

  if ( HI%isc-hshift<HI%isd .or. HI%iec+hshift>HI%ied .or. &
       HI%jscB-hshift<HI%jsdB .or. HI%jecB+hshift>HI%jedB ) then
    write(0,*) 'chksum_v_2d: haloshift =',hshift
    write(0,*) 'chksum_v_2d: isd,isc,iec,ied=',HI%isd,HI%isc,HI%iec,HI%ied
    write(0,*) 'chksum_v_2d: jsd,jsc,jec,jed=',HI%jsdB,HI%jscB,HI%jecB,HI%jedB
    call chksum_error(FATAL,'Error in chksum_v_2d '//trim(mesg))
  endif

  bc0=subchk(array, HI, 0, 0)

  if (hshift==0) then
      if (is_root_pe()) call chk_sum_msg("v-point:",bc0,mesg)
      return
  endif

  bcSW=subchk(array, HI, -hshift, -hshift)
  bcSE=subchk(array, HI, hshift, -hshift)
  bcNW=subchk(array, HI, -hshift, hshift)
  bcNE=subchk(array, HI, hshift, hshift)

  if (is_root_pe()) call chk_sum_msg("v-point:",bc0,bcSW,bcSE,bcNW,bcNE,mesg)

  contains

  integer function subchk(array, HI, di, dj)
    type(hor_index_type), intent(in) :: HI
    real, dimension(HI%isd:,HI%JsdB:), intent(in) :: array
    integer, intent(in) :: di, dj
    integer :: bitcount, i, j, bc
    subchk = 0
    do j=HI%jscB+dj,HI%jecB+dj; do i=HI%isc+di,HI%iec+di
        bc = bitcount(abs(array(i,j)))
        subchk = subchk + bc
    enddo; enddo
    call sum_across_PEs(subchk)
    subchk=mod(subchk,1000000000)
  end function subchk

  subroutine subStats(HI, array, mesg)
    type(hor_index_type), intent(in) :: HI
    real, dimension(HI%isd:,HI%JsdB:), intent(in) :: array
    character(len=*), intent(in) :: mesg
    integer :: i, j, n
    real :: aMean, aMin, aMax

    aMin = array(HI%isc,HI%jscB)
    aMax = array(HI%isc,HI%jscB)
    n = 0
    do j=HI%jscB,HI%jecB ; do i=HI%isc,HI%iec
      aMin = min(aMin, array(i,j))
      aMax = max(aMax, array(i,j))
      n = n + 1
    enddo ; enddo

    aMean = reproducing_sum(array(HI%isc:HI%iec,HI%jscB:HI%jecB))
    call sum_across_PEs(n)
    call min_across_PEs(aMin)
    call max_across_PEs(aMax)
    aMean = aMean / real(n)
    if (is_root_pe()) call chk_sum_msg("v-point:",aMean,aMin,aMax,mesg)
  end subroutine subStats

end subroutine chksum_v_2d

! =====================================================================

!> chksum_h_3d performs checksums on a 3d array staggered at tracer points.
subroutine chksum_h_3d(array, mesg, HI, haloshift, fname)
  type(hor_index_type),             intent(in) :: HI !< A horizontal index type
  real, dimension(HI%isd:,HI%jsd:,:),  intent(in) :: array !< The array to be checksummed
  character(len=*),                  intent(in) :: mesg  !< An identifying message
  integer,                 optional, intent(in) :: haloshift !< The width of halos to check (default 0)
  character(len=*),        optional, intent(in) :: fname  !< Name of file to dump

  integer :: bc0,bcSW,bcSE,bcNW,bcNE,hshift

  if (present(fname) .and. transform_input_test) then
      if (transform_input_on_this_pe) then
        call write_to_netcdf_3d(array, fname//'_h_trans.nc')
      else
        call write_to_netcdf_3d(array, fname//'_h.nc')
      endif
  endif

  if (checkForNaNs) then
    if (is_NaN(array(HI%isc:HI%iec,HI%jsc:HI%jec,:))) &
      call chksum_error(FATAL, 'NaN detected: '//trim(mesg))
!   if (is_NaN(array)) &
!     call chksum_error(FATAL, 'NaN detected in halo: '//trim(mesg))
  endif

  if (calculateStatistics) call subStats(HI, array, mesg)

  if (.not.writeChksums) return

  hshift=default_shift
  if (present(haloshift)) hshift=haloshift
  if (hshift<0) hshift=HI%ied-HI%iec

  if ( HI%isc-hshift<HI%isd .or. HI%iec+hshift>HI%ied .or. &
       HI%jsc-hshift<HI%jsd .or. HI%jec+hshift>HI%jed ) then
    write(0,*) 'chksum_h_3d: haloshift =',hshift
    write(0,*) 'chksum_h_3d: isd,isc,iec,ied=',HI%isd,HI%isc,HI%iec,HI%ied
    write(0,*) 'chksum_h_3d: jsd,jsc,jec,jed=',HI%jsd,HI%jsc,HI%jec,HI%jed
    call chksum_error(FATAL,'Error in chksum_h_3d '//trim(mesg))
  endif

  bc0=subchk(array, HI, 0, 0)

  if (hshift==0) then
      if (is_root_pe()) call chk_sum_msg("h-point:",bc0,mesg)
      return
  endif

  bcSW=subchk(array, HI, -hshift, -hshift)
  bcSE=subchk(array, HI, hshift, -hshift)
  bcNW=subchk(array, HI, -hshift, hshift)
  bcNE=subchk(array, HI, hshift, hshift)

  if (is_root_pe()) call chk_sum_msg("h-point:",bc0,bcSW,bcSE,bcNW,bcNE,mesg)

  contains

  integer function subchk(array, HI, di, dj)
    type(hor_index_type), intent(in) :: HI
    real, dimension(HI%isd:,HI%jsd:,:), intent(in) :: array
    integer, intent(in) :: di, dj
    integer :: bitcount, i, j, k, bc
    subchk = 0
    do k=LBOUND(array,3),UBOUND(array,3) ; do j=HI%jsc+dj,HI%jec+dj ; do i=HI%isc+di,HI%iec+di
      bc = bitcount(abs(array(i,j,k)))
      subchk = subchk + bc
    enddo ; enddo ; enddo
    call sum_across_PEs(subchk)
    subchk=mod(subchk,1000000000)
  end function subchk

  subroutine subStats(HI, array, mesg)
    type(hor_index_type), intent(in) :: HI
    real, dimension(HI%isd:,HI%jsd:,:), intent(in) :: array
    character(len=*), intent(in) :: mesg
    integer :: i, j, k, n
    real :: aMean, aMin, aMax

    aMin = array(HI%isc,HI%jsc,1)
    aMax = array(HI%isc,HI%jsc,1)
    n = 0
    do k=LBOUND(array,3),UBOUND(array,3) ; do j=HI%jsc,HI%jec ; do i=HI%isc,HI%iec
      aMin = min(aMin, array(i,j,k))
      aMax = max(aMax, array(i,j,k))
      n = n + 1
    enddo ; enddo ; enddo
    aMean = reproducing_sum(array(HI%isc:HI%iec,HI%jsc:HI%jec,:))
    call sum_across_PEs(n)
    call min_across_PEs(aMin)
    call max_across_PEs(aMax)
    aMean = aMean / real(n)
    if (is_root_pe()) call chk_sum_msg("h-point:",aMean,aMin,aMax,mesg)
  end subroutine subStats

end subroutine chksum_h_3d

! =====================================================================

!> chksum_B_3d performs checksums on a 3d array staggered at corner points.
subroutine chksum_B_3d(array, mesg, HI, haloshift, fname)
  type(hor_index_type),              intent(in) :: HI !< A horizontal index type
  real, dimension(HI%IsdB:,HI%JsdB:,:), intent(in) :: array !< The array to be checksummed
  character(len=*),                   intent(in) :: mesg  !< An identifying message
  integer,                  optional, intent(in) :: haloshift !< The width of halos to check (default 0)
  character(len=*),        optional, intent(in) :: fname  !< Name of file to dump

  integer :: bc0,bcSW,bcSE,bcNW,bcNE,hshift

  if (present(fname) .and. transform_input_test) then
      if (transform_input_on_this_pe) then
        call write_to_netcdf_3d(array, fname//'_B_trans.nc')
      else
        call write_to_netcdf_3d(array, fname//'_B.nc')
      endif
  endif

  if (checkForNaNs) then
    if (is_NaN(array(HI%IscB:HI%IecB,HI%JscB:HI%JecB,:))) &
      call chksum_error(FATAL, 'NaN detected: '//trim(mesg))
!   if (is_NaN(array)) &
!     call chksum_error(FATAL, 'NaN detected in halo: '//trim(mesg))
  endif

  if (calculateStatistics) call subStats(HI, array, mesg)

  if (.not.writeChksums) return

  hshift=default_shift
  if (present(haloshift)) hshift=haloshift
  if (hshift<0) hshift=HI%ied-HI%iec

  if ( HI%isc-hshift<HI%isd .or. HI%iec+hshift>HI%ied .or. &
       HI%jsc-hshift<HI%jsd .or. HI%jec+hshift>HI%jed ) then
    write(0,*) 'chksum_B_3d: haloshift =',hshift
    write(0,*) 'chksum_B_3d: isd,isc,iec,ied=',HI%isd,HI%isc,HI%iec,HI%ied
    write(0,*) 'chksum_B_3d: jsd,jsc,jec,jed=',HI%jsd,HI%jsc,HI%jec,HI%jed
    call chksum_error(FATAL,'Error in chksum_B_3d '//trim(mesg))
  endif

  bc0=subchk(array, HI, 0, 0)

  if (hshift==0) then
    if (is_root_pe()) call chk_sum_msg("B-point:",bc0,mesg)
    return
  endif

  bcSW=subchk(array, HI, -hshift, -hshift)
  bcSE=subchk(array, HI, hshift, -hshift)
  bcNW=subchk(array, HI, -hshift, hshift)
  bcNE=subchk(array, HI, hshift, hshift)

  if (is_root_pe()) call chk_sum_msg("B-point:",bc0,bcSW,bcSE,bcNW,bcNE,mesg)

  contains

  integer function subchk(array, HI, di, dj)
    type(hor_index_type), intent(in) :: HI
    real, dimension(HI%IsdB:,HI%JsdB:,:), intent(in) :: array
    integer, intent(in) :: di, dj
    integer :: bitcount, i, j, k, bc
    subchk = 0
    do k=LBOUND(array,3),UBOUND(array,3) ; do j=HI%jsc+dj,HI%jec+dj ; do i=HI%isc+di,HI%iec+di
      bc = bitcount(abs(array(i,j,k)))
      subchk = subchk + bc
    enddo ; enddo ; enddo
    call sum_across_PEs(subchk)
    subchk=mod(subchk,1000000000)
  end function subchk

  subroutine subStats(HI, array, mesg)
    type(hor_index_type), intent(in) :: HI
    real, dimension(HI%IsdB:,HI%JsdB:,:), intent(in) :: array
    character(len=*), intent(in) :: mesg
    integer :: i, j, k, n
    real :: aMean, aMin, aMax

    aMin = array(HI%isc,HI%jsc,1)
    aMax = array(HI%isc,HI%jsc,1)
    n = 0
    do k=LBOUND(array,3),UBOUND(array,3) ; do j=HI%jsc,HI%jec ; do i=HI%isc,HI%iec
      aMin = min(aMin, array(i,j,k))
      aMax = max(aMax, array(i,j,k))
      n = n + 1
    enddo ; enddo ; enddo
    aMean = reproducing_sum(array(HI%isc:HI%iec,HI%jsc:HI%jec,:))
    call sum_across_PEs(n)
    call min_across_PEs(aMin)
    call max_across_PEs(aMax)
    aMean = aMean / real(n)
    if (is_root_pe()) call chk_sum_msg("B-point:",aMean,aMin,aMax,mesg)
  end subroutine subStats

end subroutine chksum_B_3d

! =====================================================================

!> chksum_u_3d performs checksums on a 3d array staggered at C-grid u points.
subroutine chksum_u_3d(array, mesg, HI, haloshift, fname)
  type(hor_index_type),             intent(in) :: HI !< A horizontal index type
  real, dimension(HI%isdB:,HI%Jsd:,:), intent(in) :: array !< The array to be checksummed
  character(len=*),                   intent(in) :: mesg  !< An identifying message
  integer,    optional, intent(in) :: haloshift !< The width of halos to check (default 0)
  character(len=*),        optional, intent(in) :: fname  !< Name of file to dump

  integer :: bc0,bcSW,bcSE,bcNW,bcNE,hshift

  if (present(fname) .and. transform_input_test) then
      if (transform_input_on_this_pe) then
        call write_to_netcdf_3d(array, fname//'_u_trans.nc')
      else
        call write_to_netcdf_3d(array, fname//'_u.nc')
      endif
  endif

  if (checkForNaNs) then
    if (is_NaN(array(HI%IscB:HI%IecB,HI%jsc:HI%jec,:))) &
      call chksum_error(FATAL, 'NaN detected: '//trim(mesg))
!   if (is_NaN(array)) &
!     call chksum_error(FATAL, 'NaN detected in halo: '//trim(mesg))
  endif

  if (calculateStatistics) call subStats(HI, array, mesg)

  if (.not.writeChksums) return

  hshift=default_shift
  if (present(haloshift)) hshift=haloshift
  if (hshift<0) hshift=HI%ied-HI%iec

  if ( HI%isc-hshift<HI%isd .or. HI%iec+hshift>HI%ied .or. &
       HI%jsc-hshift<HI%jsd .or. HI%jec+hshift>HI%jed ) then
    write(0,*) 'chksum_u_3d: haloshift =',hshift
    write(0,*) 'chksum_u_3d: isd,isc,iec,ied=',HI%isd,HI%isc,HI%iec,HI%ied
    write(0,*) 'chksum_u_3d: jsd,jsc,jec,jed=',HI%jsd,HI%jsc,HI%jec,HI%jed
    call chksum_error(FATAL,'Error in chksum_u_3d '//trim(mesg))
  endif

  bc0=subchk(array, HI, 0, 0)

  if (hshift==0) then
    if (is_root_pe()) call chk_sum_msg("u-point:",bc0,mesg)
    return
  endif

  bcSW=subchk(array, HI, -hshift, -hshift)
  bcSE=subchk(array, HI, hshift, -hshift)
  bcNW=subchk(array, HI, -hshift, hshift)
  bcNE=subchk(array, HI, hshift, hshift)

  if (is_root_pe()) call chk_sum_msg("u-point:",bc0,bcSW,bcSE,bcNW,bcNE,mesg)

  contains

  integer function subchk(array, HI, di, dj)
    type(hor_index_type), intent(in) :: HI
    real, dimension(HI%IsdB:,HI%jsd:,:), intent(in) :: array
    integer, intent(in) :: di, dj
    integer :: bitcount, i, j, k, bc
    subchk = 0
    do k=LBOUND(array,3),UBOUND(array,3) ; do j=HI%jsc+dj,HI%jec+dj ; do i=HI%isc+di,HI%iec+di
      bc = bitcount(abs(array(i,j,k)))
      subchk = subchk + bc
    enddo ; enddo ; enddo
    call sum_across_PEs(subchk)
    subchk=mod(subchk,1000000000)
  end function subchk

  subroutine subStats(HI, array, mesg)
    type(hor_index_type), intent(in) :: HI
    real, dimension(HI%IsdB:,HI%jsd:,:), intent(in) :: array
    character(len=*), intent(in) :: mesg
    integer :: i, j, k, n
    real :: aMean, aMin, aMax

    aMin = array(HI%isc,HI%jsc,1)
    aMax = array(HI%isc,HI%jsc,1)
    n = 0
    do k=LBOUND(array,3),UBOUND(array,3) ; do j=HI%jsc,HI%jec ; do i=HI%isc,HI%iec
      aMin = min(aMin, array(i,j,k))
      aMax = max(aMax, array(i,j,k))
      n = n + 1
    enddo ; enddo ; enddo
    aMean = reproducing_sum(array(HI%isc:HI%iec,HI%jsc:HI%jec,:))
    call sum_across_PEs(n)
    call min_across_PEs(aMin)
    call max_across_PEs(aMax)
    aMean = aMean / real(n)
    if (is_root_pe()) call chk_sum_msg("u-point:",aMean,aMin,aMax,mesg)
  end subroutine subStats

end subroutine chksum_u_3d

! =====================================================================

!> chksum_v_3d performs checksums on a 3d array staggered at C-grid v points.
subroutine chksum_v_3d(array, mesg, HI, haloshift, fname)
  type(hor_index_type),             intent(in) :: HI !< A horizontal index type
  real, dimension(HI%isd:,HI%JsdB:,:), intent(in) :: array !< The array to be checksummed
  character(len=*),                  intent(in) :: mesg  !< An identifying message
  integer,    optional, intent(in) :: haloshift !< The width of halos to check (default 0)
  character(len=*),        optional, intent(in) :: fname  !< Name of file to dump

  integer :: bc0,bcSW,bcSE,bcNW,bcNE,hshift

  if (present(fname) .and. transform_input_test) then
      if (transform_input_on_this_pe) then
        call write_to_netcdf_3d(array, fname//'_v_trans.nc')
      else
        call write_to_netcdf_3d(array, fname//'_v.nc')
      endif
  endif

  if (checkForNaNs) then
    if (is_NaN(array(HI%isc:HI%iec,HI%JscB:HI%JecB,:))) &
      call chksum_error(FATAL, 'NaN detected: '//trim(mesg))
!   if (is_NaN(array)) &
!     call chksum_error(FATAL, 'NaN detected in halo: '//trim(mesg))
  endif

  if (calculateStatistics) call subStats(HI, array, mesg)

  if (.not.writeChksums) return

  hshift=default_shift
  if (present(haloshift)) hshift=haloshift
  if (hshift<0) hshift=HI%ied-HI%iec

  if ( HI%isc-hshift<HI%isd .or. HI%iec+hshift>HI%ied .or. &
       HI%jsc-hshift<HI%jsd .or. HI%jec+hshift>HI%jed ) then
    write(0,*) 'chksum_v_3d: haloshift =',hshift
    write(0,*) 'chksum_v_3d: isd,isc,iec,ied=',HI%isd,HI%isc,HI%iec,HI%ied
    write(0,*) 'chksum_v_3d: jsd,jsc,jec,jed=',HI%jsd,HI%jsc,HI%jec,HI%jed
    call chksum_error(FATAL,'Error in chksum_v_3d '//trim(mesg))
  endif

  bc0=subchk(array, HI, 0, 0)

  if (hshift==0) then
    if (is_root_pe()) call chk_sum_msg("v-point:",bc0,mesg)
    return
  endif

  bcSW=subchk(array, HI, -hshift, -hshift)
  bcSE=subchk(array, HI, hshift, -hshift)
  bcNW=subchk(array, HI, -hshift, hshift)
  bcNE=subchk(array, HI, hshift, hshift)

  if (is_root_pe()) call chk_sum_msg("v-point:",bc0,bcSW,bcSE,bcNW,bcNE,mesg)

  contains

  integer function subchk(array, HI, di, dj)
    type(hor_index_type), intent(in) :: HI
    real, dimension(HI%isd:,HI%JsdB:,:), intent(in) :: array
    integer, intent(in) :: di, dj
    integer :: bitcount, i, j, k, bc
    subchk = 0
    do k=LBOUND(array,3),UBOUND(array,3) ; do j=HI%jsc+dj,HI%jec+dj ; do i=HI%isc+di,HI%iec+di
      bc = bitcount(abs(array(i,j,k)))
      subchk = subchk + bc
    enddo ; enddo ; enddo
    call sum_across_PEs(subchk)
    subchk=mod(subchk,1000000000)
  end function subchk

  subroutine subStats(HI, array, mesg)
    type(hor_index_type), intent(in) :: HI
    real, dimension(HI%isd:,HI%JsdB:,:), intent(in) :: array
    character(len=*), intent(in) :: mesg
    integer :: i, j, k, n
    real :: aMean, aMin, aMax

    aMin = array(HI%isc,HI%jsc,1)
    aMax = array(HI%isc,HI%jsc,1)
    n = 0
    do k=LBOUND(array,3),UBOUND(array,3) ; do j=HI%jsc,HI%jec ; do i=HI%isc,HI%iec
      aMin = min(aMin, array(i,j,k))
      aMax = max(aMax, array(i,j,k))
      n = n + 1
    enddo ; enddo ; enddo
    aMean = reproducing_sum(array(HI%isc:HI%iec,HI%jsc:HI%jec,:))
    call sum_across_PEs(n)
    call min_across_PEs(aMin)
    call max_across_PEs(aMax)
    aMean = aMean / real(n)
    if (is_root_pe()) call chk_sum_msg("v-point:",aMean,aMin,aMax,mesg)
  end subroutine subStats

end subroutine chksum_v_3d


! =====================================================================

!   These are the older version of chksum that do not take the grid staggering
! into account.

!> chksum1d does a checksum of a 1-dimensional array.
subroutine chksum1d(array, mesg, start_i, end_i, compare_PEs)
  real, dimension(:), intent(in) :: array   !< The array to be summed (index starts at 1).
  character(len=*),   intent(in) :: mesg    !< An identifying message.
  integer, optional,  intent(in) :: start_i !< The starting index for the sum (default 1)
  integer, optional,  intent(in) :: end_i   !< The ending index for the sum (default all)
  logical, optional,  intent(in) :: compare_PEs !< If true, compare across PEs instead of summing
                                                !! and list the root_PE value (default true)

  integer :: is, ie, i, bc, sum1, sum_bc
  integer :: bitcount
  real :: sum
  real, allocatable :: sum_here(:)
  logical :: compare
  integer :: pe_num   ! pe number of the data
  integer :: nPEs     ! Total number of processsors

  is = LBOUND(array,1) ; ie = UBOUND(array,1)
  if (present(start_i)) is = start_i
  if (present(end_i)) ie = end_i
  compare = .true. ; if (present(compare_PEs)) compare = compare_PEs

  sum = 0.0 ; sum_bc = 0
  do i=is,ie
    sum = sum + array(i)
    bc = bitcount(ABS(array(i)))
    sum_bc = sum_bc + bc
  enddo

  pe_num = pe_here() + 1 - root_pe() ; nPEs = num_pes()
  allocate(sum_here(nPEs)) ; sum_here(:) = 0.0 ; sum_here(pe_num) = sum
  call sum_across_PEs(sum_here,nPEs)

  sum1 = sum_bc
  call sum_across_PEs(sum1)

  if (.not.compare) then
    sum = 0.0
    do i=1,nPEs ; sum = sum + sum_here(i) ; enddo
    sum_bc = sum1
  elseif (is_root_pe()) then
    if (sum1 /= nPEs*sum_bc) &
      write(0, '(A40," bitcounts do not match across PEs: ",I12,1X,I12)') &
            mesg, sum1, nPEs*sum_bc
    do i=1,nPEs ; if (sum /= sum_here(i)) then
      write(0, '(A40," PE ",i4," sum mismatches root_PE: ",3(ES22.13,1X))') &
            mesg, i, sum_here(i), sum, sum_here(i)-sum
    endif ; enddo
  endif
  deallocate(sum_here)

  if (is_root_pe()) &
    write(0,'(A50,1X,ES25.16,1X,I12)') mesg, sum, sum_bc

end subroutine chksum1d

! =====================================================================
!   These are the older version of chksum that do not take the grid staggering
! into account.

!> chksum2d does a checksum of all data in a 2-d array.
subroutine chksum2d(array, mesg)

  real, dimension(:,:) :: array
  character(len=*) :: mesg

  integer :: bitcount
  integer :: xs,xe,ys,ye,i,j,sum1,bc
  real :: sum

  xs = LBOUND(array,1) ; xe = UBOUND(array,1)
  ys = LBOUND(array,2) ; ye = UBOUND(array,2)

  sum = 0.0 ; sum1 = 0
  do i=xs,xe ; do j=ys,ye
    bc = bitcount(abs(array(i,j)))
    sum1 = sum1 + bc
  enddo ; enddo
  call sum_across_PEs(sum1)

  sum = reproducing_sum(array(:,:))

  if (is_root_pe()) &
    write(0,'(A50,1X,ES25.16,1X,I12)') mesg, sum, sum1
!    write(0,'(A40,1X,Z16.16,1X,Z16.16,1X,ES25.16,1X,I12)') &
!      mesg, sum, sum1, sum, sum1

end subroutine chksum2d

!> chksum3d does a checksum of all data in a 2-d array.
subroutine chksum3d(array, mesg)

  real, dimension(:,:,:) :: array
  character(len=*) :: mesg

  integer :: bitcount
  integer :: xs,xe,ys,ye,zs,ze,i,j,k, bc,sum1
  real :: sum

  xs = LBOUND(array,1) ; xe = UBOUND(array,1)
  ys = LBOUND(array,2) ; ye = UBOUND(array,2)
  zs = LBOUND(array,3) ; ze = UBOUND(array,3)

  sum = 0.0 ; sum1 = 0
  do i=xs,xe ; do j=ys,ye ; do k=zs,ze
    bc = bitcount(ABS(array(i,j,k)))
    sum1 = sum1 + bc
  enddo ; enddo ; enddo

  call sum_across_PEs(sum1)
  sum = reproducing_sum(array(:,:,:))

  if (is_root_pe()) &
    write(0,'(A50,1X,ES25.16,1X,I12)') mesg, sum, sum1
!    write(0,'(A40,1X,Z16.16,1X,Z16.16,1X,ES25.16,1X,I12)') &
!      mesg, sum, sum1, sum, sum1

end subroutine chksum3d

! =====================================================================

!> This function returns .true. if x is a NaN, and .false. otherwise.
function is_NaN_0d(x)
  real, intent(in) :: x !< The value to be checked for NaNs.
  logical :: is_NaN_0d

 !is_NaN_0d = (((x < 0.0) .and. (x >= 0.0)) .or. &
 !          (.not.(x < 0.0) .and. .not.(x >= 0.0)))
  if (((x < 0.0) .and. (x >= 0.0)) .or. &
            (.not.(x < 0.0) .and. .not.(x >= 0.0))) then
    is_NaN_0d = .true.
  else
    is_NaN_0d = .false.
  endif

end function is_NaN_0d

! =====================================================================

!> This function returns .true. if any element of x is a NaN, and .false. otherwise.
function is_NaN_1d(x, skip_mpp)
  real, dimension(:), intent(in) :: x !< The array to be checked for NaNs.
  logical :: is_NaN_1d
  logical, optional :: skip_mpp  !< If true, only check this array only on the local PE (default false).

  integer :: i, n
  logical :: call_mpp

  n = 0
  do i = LBOUND(x,1), UBOUND(x,1)
    if (is_NaN_0d(x(i))) n = n + 1
  enddo
  call_mpp = .true.
  if (present(skip_mpp)) call_mpp = .not.skip_mpp

  if (call_mpp) call sum_across_PEs(n)
  is_NaN_1d = .false.
  if (n>0) is_NaN_1d = .true.

end function is_NaN_1d

! =====================================================================

!> This function returns .true. if any element of x is a NaN, and .false. otherwise.
function is_NaN_2d(x)
  real, dimension(:,:), intent(in) :: x !< The array to be checked for NaNs.
  logical :: is_NaN_2d

  integer :: i, j, n

  n = 0
  do j = LBOUND(x,2), UBOUND(x,2) ; do i = LBOUND(x,1), UBOUND(x,1)
    if (is_NaN_0d(x(i,j))) n = n + 1
  enddo ; enddo
  call sum_across_PEs(n)
  is_NaN_2d = .false.
  if (n>0) is_NaN_2d = .true.

end function is_NaN_2d

! =====================================================================

!> This function returns .true. if any element of x is a NaN, and .false. otherwise.
function is_NaN_3d(x)
  real, dimension(:,:,:), intent(in) :: x !< The array to be checked for NaNs.
  logical :: is_NaN_3d

  integer :: i, j, k, n

  n = 0
  do k = LBOUND(x,3), UBOUND(x,3)
    do j = LBOUND(x,2), UBOUND(x,2) ; do i = LBOUND(x,1), UBOUND(x,1)
      if (is_NaN_0d(x(i,j,k))) n = n + 1
    enddo ; enddo
  enddo
  call sum_across_PEs(n)
  is_NaN_3d = .false.
  if (n>0) is_NaN_3d = .true.

end function is_NaN_3d

! =====================================================================

subroutine chk_sum_msg1(fmsg,bc0,mesg)
  character(len=*), intent(in) :: fmsg, mesg
  integer,          intent(in) :: bc0
  if (is_root_pe()) write(0,'(A,1(A,I10,X),A)') fmsg," c=",bc0,trim(mesg)
end subroutine chk_sum_msg1

! =====================================================================

subroutine chk_sum_msg5(fmsg,bc0,bcSW,bcSE,bcNW,bcNE,mesg)
  character(len=*), intent(in) :: fmsg, mesg
  integer,          intent(in) :: bc0,bcSW,bcSE,bcNW,bcNE
  if (is_root_pe()) write(0,'(A,5(A,I10,1X),A)') &
     fmsg," c=",bc0,"sw=",bcSW,"se=",bcSE,"nw=",bcNW,"ne=",bcNE,trim(mesg)
end subroutine chk_sum_msg5

! =====================================================================

subroutine chk_sum_msg2(fmsg,bc0,bcSW,mesg)
  character(len=*), intent(in) :: fmsg, mesg
  integer,          intent(in) :: bc0,bcSW
  if (is_root_pe()) write(0,'(A,2(A,I9,1X),A)') &
     fmsg," c=",bc0,"s/w=",bcSW,trim(mesg)
end subroutine chk_sum_msg2

! =====================================================================

subroutine chk_sum_msg3(fmsg,aMean,aMin,aMax,mesg)
  character(len=*), intent(in) :: fmsg, mesg
  real,             intent(in) :: aMean,aMin,aMax
  if (is_root_pe()) write(0,'(A,3(A,ES25.16,1X),A)') &
     fmsg," mean=",aMean,"min=",aMin,"max=",aMax,trim(mesg)
end subroutine chk_sum_msg3

! =====================================================================

!> MOM_checksums_init initializes the MOM_checksums module. As it happens, the
!! only thing that it does is to log the version of this module.
subroutine MOM_checksums_init(param_file)
  type(param_file_type),   intent(in)    :: param_file
! This include declares and sets the variable "version".
#include "version_variable.h"
  character(len=40)  :: mod = "MOM_checksums" ! This module's name.

  integer, dimension(6) :: ensemble_size
  integer, dimension(:, :), allocatable :: ensemble_pelist

  call log_version(param_file, mod, version)

  call get_param(param_file, mod, "TRANSFORM_INPUT_TEST", &
                 transform_input_test, &
                 "Whether or not to run a transformation \n"//&
                 "test. This involves transposing or rotating to all \n"//&
                 "model inputs. This is a testing feature that can be \n"//&
                 "used to help find horizontal indexing errors. \n", default=.false.)

  ! Check that we're running as an ensemble. And that each ensemble member uses 1 PE.
  if (transform_input_test) then
    ensemble_size = get_ensemble_size()
    if (ensemble_size(1) == 1) then
      call MOM_error(FATAL, &
                     "TRANSFORM_INPUT_TEST: Model not within ensemble")
    endif
    if (ensemble_size(2) /= 1) then
      call MOM_error(FATAL, &
                     "TRANSFORM_INPUT_TEST: must have 1 PE per ensemble")
    endif

    allocate(ensemble_pelist(ensemble_size(1), ensemble_size(2)))
    call get_ensemble_pelist(ensemble_pelist)

    ! For this test the root PE will be the transformed run and the other
    ! will be the vanilla run.
    if (PE_here() == ensemble_pelist(1, 1)) then
      transform_input_on_this_pe = .true.
    endif

    deallocate(ensemble_pelist)

  endif

end subroutine MOM_checksums_init

! =====================================================================

subroutine chksum_error(signal, message)
  ! Wrapper for MOM_error to help place specific break points in
  ! debuggers
  integer, intent(in) :: signal
  character(len=*), intent(in) :: message
  call MOM_error(signal, message)
end subroutine chksum_error

! =====================================================================

function do_transform_input()
    logical :: do_transform_input

    do_transform_input = transform_input_on_this_pe

end function do_transform_input

subroutine transform_input_2d(array)
  real, dimension(:,:), intent(inout) :: array !< The array to be transformed

  real, allocatable, dimension(:,:) :: tmp

  if (.not. transform_input_on_this_pe) then
    return
  endif

  if (size(array, 1) /= size(array, 2)) then
    call MOM_error(FATAL, 'transform_input_2d: transform requires a square domain.')
  endif

  allocate(tmp(size(array, 1), size(array, 2)))

  if (.true.) then
      call transpose_2d(array, tmp)
  else
      ! Try a 90 degree rotation
      call rot90_2d(array, tmp, 1)
  endif

  array(:, :) = tmp(:, :)

  deallocate(tmp)

end subroutine transform_input_2d

subroutine transform_input_3d(array)
  real, dimension(:,:,:), intent(inout) :: array !< The array to be transformed

  real, allocatable, dimension(:,:,:) :: tmp
  integer :: k

  if (.not. transform_input_on_this_pe) then
    return
  endif

  allocate(tmp(size(array, 1), size(array, 2), size(array, 3)))

  ! Try a 90 degree rotation
  if (.true.) then
      call rot90_3d(array, tmp, 1)
  else
      call transpose_3d(array, tmp)
  endif

  array(:, :, :) = tmp(:, :, :)

  deallocate(tmp)

end subroutine transform_input_3d

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
  real, dimension(:,:,:), intent(inout) :: arrayIn !< The array to be rotated
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

subroutine transpose_2d(arrayIn, arrayOut)
  real, dimension(:,:), intent(in) :: arrayIn !< Array to be transposed
  real, dimension(:,:), intent(out) :: arrayOut !< Transposed array

  arrayOut(:, :) = transpose(arrayIn(:, :))

end subroutine transpose_2d

subroutine transpose_3d(arrayIn, arrayOut)
  real, dimension(:,:,:), intent(inout) :: arrayIn !< The array to be transposed
  real, dimension(:,:,:), intent(inout) :: arrayOut !< Transposed array

  integer :: k

  do k=lbound(arrayIn, 3), ubound(arrayIn, 3)
     call transpose_2d(arrayIn(:, :, k), arrayOut(:, :, k))
  enddo

end subroutine transpose_3d

subroutine transform_and_swap_input_2d(arrayA, arrayB)
  real, intent(inout), dimension(:,:) :: arrayA, arrayB

  real, allocatable, dimension(:,:) :: tmp

  if (.not. transform_input_on_this_pe) then
    return
  endif

  if (size(arrayA, 1) /= size(arrayB, 2) .or. size(arrayA, 2) /= size(arrayB, 1)) then
    call MOM_error(FATAL, 'transform_and_swap_input_2d: arrays shapes imcompatible.')
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

end subroutine transform_and_swap_input_2d

subroutine transform_and_swap_input_3d(arrayA, arrayB)
  real, intent(inout), dimension(:,:,:) :: arrayA, arrayB

  real, allocatable, dimension(:,:,:) :: tmp

  if (.not. transform_input_on_this_pe) then
    return
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

end subroutine transform_and_swap_input_3d

subroutine compare_within_ensemble(sbuf)

  real, intent(in), dimension(:) :: sbuf

  integer, dimension(:, :), allocatable :: ensemble_pelist
  integer, dimension(6) :: ensemble_size
  real, dimension(:), allocatable :: rbuf
  integer :: sbuf_size, e, i

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

  ! Gather between the root of every ensemble member.
  call mpp_gather(sbuf, rbuf, ensemble_pelist(:, 1))

  ! Check that sbuf is the same on all pes.
  if (transform_input_on_this_pe) then
    do e=0,ensemble_size(1)-1
      do i=1,sbuf_size
        if (sbuf(i) /= rbuf(e*sbuf_size + i)) then
          call MOM_error(FATAL, &
                         'TRANSFORM_INPUT_TEST failed')
        endif
      enddo
    enddo
    print*, 'TRANSFORM_INPUT_TEST passed'
  endif

  deallocate(rbuf)


end subroutine compare_within_ensemble

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

end module MOM_checksums
