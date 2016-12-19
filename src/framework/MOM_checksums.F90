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

implicit none ; private

public :: hchksum, Bchksum, uchksum, vchksum, qchksum, chksum, is_NaN
public :: write_to_netcdf, swap_md, sym_trans, sym_trans_active, uvchksum
public :: MOM_checksums_init

interface hchksum
  module procedure chksum_h_2d, chksum_h_3d
end interface

interface Bchksum
  module procedure chksum_B_2d, chksum_B_3d
end interface

interface uchksum
  module procedure chksum_u_2d, chksum_u_3d
end interface

interface vchksum
  module procedure chksum_v_2d, chksum_v_3d
end interface

interface uvchksum
  module procedure chksum_uv_2d, chksum_uv_3d
end interface

! This is an older interface that has been renamed Bchksum
interface qchksum
  module procedure chksum_B_2d, chksum_B_3d
end interface

interface chksum
  module procedure chksum1d, chksum2d, chksum3d
end interface

interface chk_sum_msg
  module procedure chk_sum_msg1, chk_sum_msg2, chk_sum_msg3, chk_sum_msg4, chk_sum_msg5
end interface

interface is_NaN
  module procedure is_NaN_0d, is_NaN_1d, is_NaN_2d, is_NaN_3d
end interface

interface swap_md
  module procedure swap_2d, swap_3d
end interface

interface sym_trans
  module procedure sym_trans_2d, sym_trans_3d
end interface

interface write_to_netcdf
  module procedure write_to_netcdf2d, write_to_netcdf3d
end interface


integer, parameter :: default_shift=0
logical :: calculateStatistics=.true. ! If true, report min, max and mean.
logical :: writeChksums=.true. ! If true, report the bitcount checksum
logical :: checkForNaNs=.true. ! If true, checks array for NaNs and cause
                               ! FATAL error is any are found
logical :: sym_trans_is_configured = .false.
logical :: sym_trans_is_active = .false.

contains

! =====================================================================

!> chksum_h_2d performs checksums on a 2d array staggered at tracer points.
subroutine chksum_h_2d(array, mesg, HI, haloshift)
  type(hor_index_type),           intent(in) :: HI     !< A horizontal index type
  real, dimension(HI%isd:,HI%jsd:), intent(in) :: array !< The array to be checksummed
  character(len=*),                intent(in) :: mesg  !< An identifying message
  integer,               optional, intent(in) :: haloshift !< The width of halos to check (default 0)

  integer :: bc0,bcSW,bcSE,bcNW,bcNE,hshift

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

  if (sym_trans_active()) then
    if (is_root_pe()) call chk_sum_msg("h-point:",bc0,bcSW,bcNW,bcSE,bcNE,mesg)
  else
    if (is_root_pe()) call chk_sum_msg("h-point:",bc0,bcSW,bcSE,bcNW,bcNE,mesg)
  endif

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

!> chksum_B_2d performs checksums on a 2d array staggered at corner points.
subroutine chksum_B_2d(array, mesg, HI, haloshift, symmetric)
  type(hor_index_type), intent(in) :: HI     !< A horizontal index type
  real, dimension(HI%IsdB:,HI%JsdB:), &
                        intent(in) :: array !< The array to be checksummed
  character(len=*),     intent(in) :: mesg  !< An identifying message
  integer,    optional, intent(in) :: haloshift !< The width of halos to check (default 0)
  logical,    optional, intent(in) :: symmetric !< If true, do the checksums on the
                                                !! full symmetric computational domain.

  integer :: bc0,bcSW,bcSE,bcNW,bcNE,hshift
  logical :: sym

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

  if (sym_trans_active()) then
    if (is_root_pe()) call chk_sum_msg("B-point:",bc0,bcSW,bcNW,bcSE,bcNE,mesg)
  else
    if (is_root_pe()) call chk_sum_msg("B-point:",bc0,bcSW,bcSE,bcNW,bcNE,mesg)
  endif

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

!> chksum_uv_2d calls checksum subroutines for both u and v points.
subroutine chksum_uv_2d(u_array, v_array, u_mesg, v_mesg, &
                        HI, uhaloshift, vhaloshift)
  real, dimension(:,:), intent(in) :: u_array !< The array to be checksummed
  real, dimension(:,:), intent(in) :: v_array !< The array to be checksummed
  character(len=*),                intent(in) :: u_mesg  !< An identifying message
  character(len=*),                intent(in) :: v_mesg  !< An identifying message
  type(hor_index_type),           intent(in) :: HI     !< A horizontal index type
  integer,               optional, intent(in) :: uhaloshift !< The width of halos to check (default 0)
  integer,               optional, intent(in) :: vhaloshift !< The width of halos to check (default 0)


  if (sym_trans_active()) then
    if (present(uhaloshift)) then
      call chksum_u_2d(v_array, u_mesg, HI, uhaloshift)
    else
      call chksum_u_2d(v_array, u_mesg, HI)
    endif
    if (present(vhaloshift)) then
      call chksum_v_2d(u_array, v_mesg, HI, vhaloshift)
    else
      call chksum_v_2d(u_array, v_mesg, HI)
    endif
  else
    if (present(uhaloshift)) then
      call chksum_u_2d(u_array, u_mesg, HI, uhaloshift)
    else
      call chksum_u_2d(u_array, u_mesg, HI)
    endif
    if (present(vhaloshift)) then
      call chksum_v_2d(v_array, v_mesg, HI, vhaloshift)
    else
      call chksum_v_2d(v_array, v_mesg, HI)
    endif
  endif

end subroutine chksum_uv_2d

!> chksum_uv_3d calls checksum subroutines for both u and v points.
subroutine chksum_uv_3d(u_array, v_array, u_mesg, v_mesg, &
                        HI, uhaloshift, vhaloshift)
  real, dimension(:,:,:), intent(in) :: u_array !< The array to be checksummed
  real, dimension(:,:,:), intent(in) :: v_array !< The array to be checksummed
  character(len=*),                intent(in) :: u_mesg  !< An identifying message
  character(len=*),                intent(in) :: v_mesg  !< An identifying message
  integer,               optional, intent(in) :: uhaloshift !< The width of halos to check (default 0)
  integer,               optional, intent(in) :: vhaloshift !< The width of halos to check (default 0)
  type(hor_index_type),           intent(in) :: HI     !< A horizontal index type


  if (sym_trans_active()) then
    if (present(uhaloshift)) then
      call chksum_u_3d(v_array, u_mesg, HI, uhaloshift)
    else
      call chksum_u_3d(v_array, u_mesg, HI)
    endif
    if (present(vhaloshift)) then
      call chksum_v_3d(u_array, v_mesg, HI, vhaloshift)
    else
      call chksum_v_3d(u_array, v_mesg, HI)
    endif
  else
    if (present(uhaloshift)) then
      call chksum_u_3d(u_array, u_mesg, HI, uhaloshift)
    else
      call chksum_u_3d(u_array, u_mesg, HI)
    endif
    if (present(vhaloshift)) then
      call chksum_v_3d(v_array, v_mesg, HI, vhaloshift)
    else
      call chksum_v_3d(v_array, v_mesg, HI)
    endif
  endif

end subroutine chksum_uv_3d

!> chksum_u_2d performs checksums on a 2d array staggered at C-grid u points.
subroutine chksum_u_2d(array, mesg, HI, haloshift)
  type(hor_index_type),           intent(in) :: HI     !< A horizontal index type
  real, dimension(HI%IsdB:,HI%jsd:), intent(in) :: array !< The array to be checksummed
  character(len=*),                intent(in) :: mesg  !< An identifying message
  integer,               optional, intent(in) :: haloshift !< The width of halos to check (default 0)

  integer :: bc0,bcSW,bcSE,bcNW,bcNE,hshift

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

  if (sym_trans_active()) then
    if (is_root_pe()) call chk_sum_msg("u-point:",bc0,bcSW,bcNW,bcSE,bcNE,mesg)
  else
    if (is_root_pe()) call chk_sum_msg("u-point:",bc0,bcSW,bcSE,bcNW,bcNE,mesg)
  endif

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
subroutine chksum_v_2d(array, mesg, HI, haloshift)
  type(hor_index_type),           intent(in) :: HI     !< A horizontal index type
  real, dimension(HI%isd:,HI%JsdB:), intent(in) :: array !< The array to be checksummed
  character(len=*),                intent(in) :: mesg  !< An identifying message
  integer,               optional, intent(in) :: haloshift !< The width of halos to check (default 0)

  integer :: bc0,bcSW,bcSE,bcNW,bcNE,hshift

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

  if (sym_trans_active()) then
    if (is_root_pe()) call chk_sum_msg("v-point:",bc0,bcSW,bcNW,bcSE,bcNE,mesg)
  else
    if (is_root_pe()) call chk_sum_msg("v-point:",bc0,bcSW,bcSE,bcNW,bcNE,mesg)
  endif

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
subroutine chksum_h_3d(array, mesg, HI, haloshift)
  type(hor_index_type),             intent(in) :: HI !< A horizontal index type
  real, dimension(HI%isd:,HI%jsd:,:),  intent(in) :: array !< The array to be checksummed
  character(len=*),                  intent(in) :: mesg  !< An identifying message
  integer,                 optional, intent(in) :: haloshift !< The width of halos to check (default 0)

  integer :: bc0,bcSW,bcSE,bcNW,bcNE,hshift

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

  if (sym_trans_active()) then
    if (is_root_pe()) call chk_sum_msg("h-point:",bc0,bcSW,bcNW,bcSE,bcNE,mesg)
  else
    if (is_root_pe()) call chk_sum_msg("h-point:",bc0,bcSW,bcSE,bcNW,bcNE,mesg)
  endif

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
    real :: aMean, aMin, aMax, hash

    hash = 0.0
    n = 1
    if (sym_trans_active()) then
      do k=LBOUND(array,3),UBOUND(array,3) ; do j=HI%jsc,HI%jec ; do i=HI%isc,HI%iec
        hash = hash + array(j,i,k)*n
        n = n + 1
      enddo ; enddo ; enddo
    else
      do k=LBOUND(array,3),UBOUND(array,3) ; do j=HI%jsc,HI%jec ; do i=HI%isc,HI%iec
        hash = hash + array(i,j,k)*n
        n = n + 1
      enddo ; enddo ; enddo
    endif

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
    if (is_root_pe()) call chk_sum_msg("h-point:",aMean,aMin,aMax,hash,mesg)
  end subroutine subStats

end subroutine chksum_h_3d

! =====================================================================

!> chksum_B_3d performs checksums on a 3d array staggered at corner points.
subroutine chksum_B_3d(array, mesg, HI, haloshift)
  type(hor_index_type),              intent(in) :: HI !< A horizontal index type
  real, dimension(HI%IsdB:,HI%JsdB:,:), intent(in) :: array !< The array to be checksummed
  character(len=*),                   intent(in) :: mesg  !< An identifying message
  integer,                  optional, intent(in) :: haloshift !< The width of halos to check (default 0)

  integer :: bc0,bcSW,bcSE,bcNW,bcNE,hshift

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

  if (sym_trans_active()) then
    if (is_root_pe()) call chk_sum_msg("B-point:",bc0,bcSW,bcNW,bcSE,bcNE,mesg)
  else
    if (is_root_pe()) call chk_sum_msg("B-point:",bc0,bcSW,bcSE,bcNW,bcNE,mesg)
  endif

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
subroutine chksum_u_3d(array, mesg, HI, haloshift)
  type(hor_index_type),             intent(in) :: HI !< A horizontal index type
  real, dimension(HI%isdB:,HI%Jsd:,:), intent(in) :: array !< The array to be checksummed
  character(len=*),                  intent(in) :: mesg  !< An identifying message
  integer,    optional, intent(in) :: haloshift !< The width of halos to check (default 0)

  integer :: bc0,bcSW,bcSE,bcNW,bcNE,hshift

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

  if (sym_trans_active()) then
    if (is_root_pe()) call chk_sum_msg("u-point:",bc0,bcSW,bcNW,bcSE,bcNE,mesg)
  else
    if (is_root_pe()) call chk_sum_msg("u-point:",bc0,bcSW,bcSE,bcNW,bcNE,mesg)
  endif

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
subroutine chksum_v_3d(array, mesg, HI, haloshift)
  type(hor_index_type),             intent(in) :: HI !< A horizontal index type
  real, dimension(HI%isd:,HI%JsdB:,:), intent(in) :: array !< The array to be checksummed
  character(len=*),                  intent(in) :: mesg  !< An identifying message
  integer,    optional, intent(in) :: haloshift !< The width of halos to check (default 0)

  integer :: bc0,bcSW,bcSE,bcNW,bcNE,hshift

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

  if (sym_trans_active()) then
    if (is_root_pe()) call chk_sum_msg("v-point:",bc0,bcSW,bcNW,bcSE,bcNE,mesg)
  else
    if (is_root_pe()) call chk_sum_msg("v-point:",bc0,bcSW,bcSE,bcNW,bcNE,mesg)
  endif

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

subroutine chksum2d(array, mesg, start_x, end_x, start_y, end_y)

  real, dimension(:,:) :: array
  character(len=*) :: mesg
  integer, optional :: start_x, end_x, start_y, end_y

  integer :: bitcount
  integer :: xs,xe,ys,ye,i,j,sum1,bc
  real :: sum

  xs = LBOUND(array,1) ; xe = UBOUND(array,1)
  ys = LBOUND(array,2) ; ye = UBOUND(array,2)
  if (present(start_x)) xs = start_x
  if (present(end_x  )) xe = end_x
  if (present(start_y)) ys = start_y
  if (present(end_y  )) ye = end_y

  sum = 0.0 ; sum1 = 0
  do i=xs,xe ; do j=ys,ye
    bc = bitcount(abs(array(i,j)))
    sum1 = sum1 + bc
  enddo ; enddo
  call sum_across_PEs(sum1)

  sum = reproducing_sum(array(xs:xe,ys:ye))

  if (is_root_pe()) &
    write(0,'(A50,1X,ES25.16,1X,I12)') mesg, sum, sum1
!    write(0,'(A40,1X,Z16.16,1X,Z16.16,1X,ES25.16,1X,I12)') &
!      mesg, sum, sum1, sum, sum1

end subroutine chksum2d

subroutine chksum3d(array, mesg, start_x, end_x, start_y, end_y, start_z, end_z)

  real, dimension(:,:,:) :: array
  character(len=*) :: mesg
  integer, optional :: start_x, end_x, start_y, end_y, start_z, end_z

  integer :: bitcount
  integer :: xs,xe,ys,ye,zs,ze,i,j,k, bc,sum1
  real :: sum

  xs = LBOUND(array,1) ; xe = UBOUND(array,1)
  ys = LBOUND(array,2) ; ye = UBOUND(array,2)
  zs = LBOUND(array,3) ; ze = UBOUND(array,3)
  if (present(start_x)) xs = start_x
  if (present(end_x  )) xe = end_x
  if (present(start_y)) ys = start_y
  if (present(end_y  )) ye = end_y
  if (present(start_z)) zs = start_z
  if (present(end_z  )) ze = end_z

  sum = 0.0 ; sum1 = 0
  do i=xs,xe ; do j=ys,ye ; do k=zs,ze
    bc = bitcount(ABS(array(i,j,k)))
    sum1 = sum1 + bc
  enddo ; enddo ; enddo

  call sum_across_PEs(sum1)
  sum = reproducing_sum(array(xs:xe,ys:ye,zs:ze))

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

subroutine chk_sum_msg4(fmsg,aMean,aMin,aMax,hash,mesg)
  character(len=*), intent(in) :: fmsg, mesg
  real,             intent(in) :: aMean,aMin,aMax,hash
  if (is_root_pe()) write(0,'(A,4(A,ES25.16,1X),A)') &
     fmsg," mean=",aMean,"min=",aMin,"max=",aMax,"hash=",hash,trim(mesg)
end subroutine chk_sum_msg4

! =====================================================================

!> MOM_checksums_init initializes the MOM_checksums module. As it happens, the
!! only thing that it does is to log the version of this module.
subroutine MOM_checksums_init(param_file)
  type(param_file_type),   intent(in)    :: param_file
! This include declares and sets the variable "version".
#include "version_variable.h"
  character(len=40)  :: mod = "MOM_checksums" ! This module's name.

  call log_version(param_file, mod, version)

  call get_param(param_file, mod, "APPLY_SYMMETRIC_INPUT_TRANSFORM", &
                 sym_trans_is_configured, &
                 "Whether or not to apply a symmetric transformation to all \n"//&
                 "model inputs. This is a testing feature that can be \n"//&
                 "used to help find horizontal indexing errors. \n", default=.false.)
  sym_trans_is_active = .false.

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

function sym_trans_active()
    logical :: sym_trans_active

    sym_trans_active = sym_trans_is_configured

end function sym_trans_active

subroutine sym_trans_2d(array)
  real, dimension(:,:), intent(inout) :: array !< The array to be transformed

  if (.not. sym_trans_is_configured) then
    return
  endif

  sym_trans_is_active = .true.

  if (size(array, 1) /= size(array, 2)) then
    call MOM_error(FATAL, 'sym_trans_2d: transform requires a square domain.')
  endif

  ! Let's try straight transpose
  ! array = transpose(array)

  ! Try a 90 degree rotation
  call rot90_2d(array, 1)

end subroutine sym_trans_2d

subroutine sym_trans_3d(array)
  real, dimension(:,:,:), intent(inout) :: array !< The array to be transformed

  integer :: k

  if (.not. sym_trans_is_configured) then
    return
  endif

  sym_trans_is_active = .true.

  if (size(array, 1) /= size(array, 2)) then
    call MOM_error(FATAL, 'sym_trans_2d: transform requires a square domain.')
  endif

  ! Try a 90 degree rotation
  call rot90_3d(array, 1)

  !do k=lbound(array, 3), ubound(array, 3)
  !  array(:, :, k) = transpose(array(:, :, k))
  !  array(:, :, k) = array(:, ubound(array, 2):lbound(array, 2):-1, k)
  !enddo

end subroutine sym_trans_3d

subroutine rot90_2d(array, nrot90)
  real, dimension(:,:), intent(inout) :: array !< The array to be rotated
  integer, intent(in) :: nrot90 !< Number of 90 degree rotations to perform

  if (.not. nrot90 < 4) then
    call MOM_error(FATAL, 'rot90_2d: nrot should be < 4')
  endif

  if (modulo(nrot90, 2) == 1) then
    if (size(array, 1) /= size(array, 2)) then
      call MOM_error(FATAL, 'rot90_2d: 90 deg rotation requires a square domain.')
    endif
  endif

  if (nrot90 == 1) then
    ! transpose, reverse rows
    array = transpose(array)
    array(:, :) = array(:, ubound(array, 2):lbound(array, 2):-1)
  elseif (nrot90 == 2) then
    ! reverse both rows and cols
    array(:, :) = array(ubound(array, 1):lbound(array, 1):-1, &
                        ubound(array, 2):lbound(array, 2):-1)
  elseif (nrot90 == 3) then
    ! transpose, reverse cols
    array = transpose(array)
    array(:, :) = array(ubound(array, 1):lbound(array, 1):-1, :)
  endif

end subroutine rot90_2d

subroutine rot90_3d(array, nrot90)
  real, dimension(:,:,:), intent(inout) :: array !< The array to be rotated
  integer, intent(in) :: nrot90 !< Number of 90 degree rotations to perform

  integer :: k

  if (.not. nrot90 < 4) then
    call MOM_error(FATAL, 'rot90_2d: nrot should be < 4')
  endif

  if (modulo(nrot90, 2) == 1) then
    if (size(array, 1) /= size(array, 2)) then
      call MOM_error(FATAL, 'rot90_2d: 90 deg rotation requires a square domain.')
    endif
  endif

  if (nrot90 == 1) then
    ! transpose, reverse rows
    do k=lbound(array, 3), ubound(array, 3)
      array(:, :, k) = transpose(array(:, :, k))
      array(:, :, k) = array(:, ubound(array, 2):lbound(array, 2):-1, k)
    enddo
  elseif (nrot90 == 2) then
    ! reverse both rows and cols
    do k=lbound(array, 3), ubound(array, 3)
        array(:, :, k) = array(ubound(array, 1):lbound(array, 1):-1, &
                               ubound(array, 2):lbound(array, 2):-1, k)
    enddo
  elseif (nrot90 == 3) then
    ! transpose, reverse cols
    do k=lbound(array, 3), ubound(array, 3)
      array(:, :, k) = transpose(array(:, :, k))
      array(:, :, k) = array(ubound(array, 1):lbound(array, 1):-1, :, k)
    enddo
  endif

end subroutine rot90_3d

subroutine swap_2d(arrayA, arrayB)
  real, intent(inout), dimension(:,:) :: arrayA, arrayB

  real, allocatable, dimension(:,:) :: tmp

  if (.not. sym_trans_is_configured) then
    return
  endif

  if (size(arrayA, 1) /= size(arrayB, 1) .or. size(arrayA, 2) /= size(arrayB, 2)) then
    call MOM_error(FATAL, 'swap_2d: different shaped arrays cannot be swapped.')
  endif

  allocate(tmp(size(arrayA, 1), size(arrayA, 2)))

  tmp(:, :) = arrayA(:, :)
  arrayA(:, :) = arrayB(:, :)
  arrayB(:, :) = tmp(:, :)

  deallocate(tmp)

end subroutine swap_2d

subroutine swap_3d(arrayA, arrayB)
  real, intent(inout), dimension(:,:,:) :: arrayA, arrayB

  real, allocatable, dimension(:,:,:) :: tmp

  if (.not. sym_trans_is_configured) then
    return
  endif

  if (size(arrayA, 1) /= size(arrayB, 1) .or. &
          size(arrayA, 2) /= size(arrayB, 2) .or. &
          size(arrayA, 3) /= size(arrayB, 3)) then
    call MOM_error(FATAL, 'swap_3d: different shaped arrays cannot be swapped.')
  endif

  allocate(tmp(size(arrayA, 1), size(arrayA, 2), size(arrayA, 3)))

  tmp(:, :, :) = arrayA(:, :, :)
  arrayA(:, :, :) = arrayB(:, :, :)
  arrayB(:, :, :) = tmp(:, :, :)

  deallocate(tmp)

end subroutine swap_3d

subroutine write_to_netcdf3d(array, file_name)
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
  ierr = nf90_def_var(file_id, 'Array',  NF90_REAL, arrdims, array_id)

  ! ...and assign units to them as an attribute 
  ierr = nf90_put_att(file_id, array_id, "units", arrunit)

  ! done defining
  ierr = nf90_enddef(file_id)

  ! Write out the values
  ierr = nf90_put_var(file_id, array_id, array)

  ! close; done
  ierr = nf90_close(file_id)
end subroutine write_to_netcdf3d

subroutine write_to_netcdf2d(array, file_name)
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
  ierr = nf90_def_var(file_id, 'Array',  NF90_REAL, arrdims, array_id)

  ! ...and assign units to them as an attribute 
  ierr = nf90_put_att(file_id, array_id, "units", arrunit)

  ! done defining
  ierr = nf90_enddef(file_id)

  ! Write out the values
  ierr = nf90_put_var(file_id, array_id, array)

  ! close; done
  ierr = nf90_close(file_id)
end subroutine write_to_netcdf2d

end module MOM_checksums
