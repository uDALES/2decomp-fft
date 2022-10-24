!=======================================================================
! This is part of the 2DECOMP&FFT library
!
! 2DECOMP&FFT is a software framework for general-purpose 2D (pencil)
! decomposition. It also implements a highly scalable distributed
! three-dimensional Fast Fourier Transform (FFT).
!
! Copyright (C) 2009-2021 Ning Li, the Numerical Algorithms Group (NAG)
!
!=======================================================================

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! Halo cell support for neighbouring pencils to exchange data
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine update_halo_real_short(in, out, level, opt_global, opt_pencil)

    implicit none

    integer, intent(IN) :: level      ! levels of halo cells required
    real(mytype), dimension(:,:,:), intent(IN) :: in
    real(mytype), allocatable, dimension(:,:,:), intent(OUT) :: out
    logical, optional :: opt_global
    integer, intent(in), optional :: opt_pencil

    call update_halo(in, out, level, decomp_main, opt_global, opt_pencil)

  end subroutine update_halo_real_short

  subroutine update_halo_real(in, out, level, decomp, opt_global, opt_pencil)

    implicit none

    integer, intent(IN) :: level      ! levels of halo cells required
    real(mytype), dimension(:,:,:), intent(IN) :: in
    real(mytype), allocatable, dimension(:,:,:), intent(OUT) :: out
    TYPE(DECOMP_INFO), intent(in) :: decomp
    logical, optional :: opt_global
    integer, intent(in), optional :: opt_pencil

    logical :: global

    ! starting/ending index of array with halo cells
    integer :: xs, ys, zs, xe, ye, ze

    integer :: i, j, k, s1, s2, s3
    integer :: data_type

    integer :: ipencil
    logical, save :: first_call_x = .true., first_call_y = .true., first_call_z = .true.
    
    data_type = real_type

#include "halo_common.f90"

    return
  end subroutine update_halo_real

  subroutine update_halo_complex_short(in, out, level, opt_global, opt_pencil)

    implicit none

    integer, intent(IN) :: level      ! levels of halo cells required
    complex(mytype), dimension(:,:,:), intent(IN) :: in
    complex(mytype), allocatable, dimension(:,:,:), intent(OUT) :: out
    logical, optional :: opt_global
    integer, intent(in), optional :: opt_pencil

    call update_halo(in, out, level, decomp_main, opt_global, opt_pencil)

  end subroutine update_halo_complex_short

  subroutine update_halo_complex(in, out, level, decomp, opt_global, opt_pencil)

    implicit none

    integer, intent(IN) :: level      ! levels of halo cells required
    complex(mytype), dimension(:,:,:), intent(IN) :: in
    complex(mytype), allocatable, dimension(:,:,:), intent(OUT) :: out
    TYPE(DECOMP_INFO), intent(in) :: decomp
    logical, optional :: opt_global
    integer, intent(in), optional :: opt_pencil

    logical :: global

    ! starting/ending index of array with halo cells
    integer :: xs, ys, zs, xe, ye, ze

    integer :: i, j, k, s1, s2, s3
    integer :: data_type

    integer :: ipencil
    logical, save :: first_call_x = .true., first_call_y = .true., first_call_z = .true.
    
    data_type = complex_type

#include "halo_common.f90"

    return
  end subroutine update_halo_complex

  subroutine exchange_halo_x_real(inout, opt_decomp, opt_xlevel)

    implicit none

    real(mytype), dimension(:,:,:), intent(INOUT) :: inout
    TYPE(DECOMP_INFO), optional :: opt_decomp
    integer, dimension(3), optional :: opt_xlevel

    TYPE(DECOMP_INFO) :: decomp
    integer :: level_x, level_y, level_z
    integer :: ierror
    integer :: icount, ilength, ijump
    integer :: halo12
    integer, dimension(4) :: requests
    integer, dimension(MPI_STATUS_SIZE,4) :: status
    integer :: tag_n, tag_s, tag_t, tag_b
    integer :: data_type
    integer :: xs, ys, zs, ye, ze, s1, s2, s3

    data_type = real_type

#include "halo_comm_x.f90"

    return
  end subroutine exchange_halo_x_real


  subroutine exchange_halo_x_complex(inout, opt_decomp, opt_xlevel)

    implicit none

    complex(mytype), dimension(:,:,:), intent(INOUT) :: inout
    TYPE(DECOMP_INFO), optional :: opt_decomp
    integer, dimension(3), optional :: opt_xlevel

    TYPE(DECOMP_INFO) :: decomp
    integer :: level_x, level_y, level_z
    integer :: ierror
    integer :: icount, ilength, ijump
    integer :: halo12
    integer, dimension(4) :: requests
    integer, dimension(MPI_STATUS_SIZE,4) :: status
    integer :: tag_n, tag_s, tag_t, tag_b
    integer :: data_type
    integer :: xs, ys, zs, ye, ze, s1, s2, s3

    data_type = complex_type

#include "halo_comm_x.f90"

    return
  end subroutine exchange_halo_x_complex


  subroutine exchange_halo_y_real(inout, opt_decomp, opt_ylevel)

    implicit none

    real(mytype), dimension(:,:,:), intent(INOUT) :: inout
    TYPE(DECOMP_INFO), optional :: opt_decomp
    integer, dimension(3), optional :: opt_ylevel

    TYPE(DECOMP_INFO) :: decomp
    integer :: level_x, level_y, level_z
    integer :: ierror
    integer :: icount, ilength, ijump
    integer :: halo21
    integer, dimension(4) :: requests
    integer, dimension(MPI_STATUS_SIZE,4) :: status
    integer :: tag_e, tag_w, tag_t, tag_b
    integer :: data_type
    integer :: xs, ys, zs, xe, ze, s1, s2, s3

    data_type = real_type

#include "halo_comm_y.f90"

    return
  end subroutine exchange_halo_y_real


  subroutine exchange_halo_y_complex(inout, opt_decomp, opt_ylevel)

    implicit none

    complex(mytype), dimension(:,:,:), intent(INOUT) :: inout
    TYPE(DECOMP_INFO), optional :: opt_decomp
    integer, dimension(3), optional :: opt_ylevel

    TYPE(DECOMP_INFO) :: decomp
    integer :: level_x, level_y, level_z
    integer :: ierror
    integer :: icount, ilength, ijump
    integer :: halo21
    integer, dimension(4) :: requests
    integer, dimension(MPI_STATUS_SIZE,4) :: status
    integer :: tag_e, tag_w, tag_t, tag_b
    integer :: data_type
    integer :: xs, ys, zs, xe, ze, s1, s2, s3

    data_type = complex_type

#include "halo_comm_y.f90"

    return
  end subroutine exchange_halo_y_complex


  subroutine exchange_halo_z_real(inout, opt_decomp, opt_zlevel)
    implicit none
    real(mytype), dimension(:,:,:), intent(INOUT) :: inout
    TYPE(DECOMP_INFO), optional :: opt_decomp
    integer, dimension(3), optional :: opt_zlevel

    TYPE(DECOMP_INFO) :: decomp
    integer :: level_x, level_y, level_z
    integer :: ierror
    integer :: icount, ilength, ijump
    integer :: halo31, halo32
    integer, dimension(4) :: requests
    integer, dimension(MPI_STATUS_SIZE,4) :: status
    integer :: tag_e, tag_w, tag_n, tag_s
    integer :: data_type
    integer :: xs, ys, zs, xe, ye, s1, s2, s3

    data_type = real_type

#include "halo_comm_z.f90"

    return
  end subroutine exchange_halo_z_real

  subroutine exchange_halo_z_complex(inout, opt_decomp, opt_zlevel)

    implicit none

    complex(mytype), dimension(:,:,:), intent(INOUT) :: inout
    TYPE(DECOMP_INFO), optional :: opt_decomp
    integer, dimension(3), optional :: opt_zlevel

    TYPE(DECOMP_INFO) :: decomp
    integer :: level_x, level_y, level_z
    integer :: ierror
    integer :: icount, ilength, ijump
    integer :: halo31, halo32
    integer, dimension(4) :: requests
    integer, dimension(MPI_STATUS_SIZE,4) :: status
    integer :: tag_e, tag_w, tag_n, tag_s
    integer :: data_type
    integer :: xs, ys, zs, xe, ye, s1, s2, s3

    data_type = complex_type

#include "halo_comm_z.f90"

    return
  end subroutine exchange_halo_z_complex

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! To support halo-cell exchange:
  !   find the MPI ranks of neighbouring pencils
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine init_neighbour

    integer :: ierror

    ! For X-pencil
    neighbour(1,1) = MPI_PROC_NULL               ! east
    neighbour(1,2) = MPI_PROC_NULL               ! west
    call MPI_CART_SHIFT(DECOMP_2D_COMM_CART_X, 0, 1, &
         neighbour(1,4), neighbour(1,3), ierror) ! north & south
    if (ierror /= 0) call decomp_2d_abort(__FILE__, __LINE__, ierror, "MPI_CART_SHIFT")
    call MPI_CART_SHIFT(DECOMP_2D_COMM_CART_X, 1, 1, &
         neighbour(1,6), neighbour(1,5), ierror) ! top & bottom
    if (ierror /= 0) call decomp_2d_abort(__FILE__, __LINE__, ierror, "MPI_CART_SHIFT")

    ! For Y-pencil
    call MPI_CART_SHIFT(DECOMP_2D_COMM_CART_Y, 0, 1, &
         neighbour(2,2), neighbour(2,1), ierror) ! east & west
    if (ierror /= 0) call decomp_2d_abort(__FILE__, __LINE__, ierror, "MPI_CART_SHIFT")
    neighbour(2,3) = MPI_PROC_NULL               ! north
    neighbour(2,4) = MPI_PROC_NULL               ! south
    call MPI_CART_SHIFT(DECOMP_2D_COMM_CART_Y, 1, 1, &
         neighbour(2,6), neighbour(2,5), ierror) ! top & bottom
    if (ierror /= 0) call decomp_2d_abort(__FILE__, __LINE__, ierror, "MPI_CART_SHIFT")

    ! For Z-pencil
    call MPI_CART_SHIFT(DECOMP_2D_COMM_CART_Z, 0, 1, &
         neighbour(3,2), neighbour(3,1), ierror) ! east & west
    if (ierror /= 0) call decomp_2d_abort(__FILE__, __LINE__, ierror, "MPI_CART_SHIFT")
    call MPI_CART_SHIFT(DECOMP_2D_COMM_CART_Z, 1, 1, &
         neighbour(3,4), neighbour(3,3), ierror) ! north & south
    if (ierror /= 0) call decomp_2d_abort(__FILE__, __LINE__, ierror, "MPI_CART_SHIFT")
    neighbour(3,5) = MPI_PROC_NULL               ! top
    neighbour(3,6) = MPI_PROC_NULL               ! bottom
    return
  end subroutine init_neighbour
