!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! This example calculates the divergence of a random field using
!   (1) global transposition
!   (2) halo-cell exchange
! The two method should give identical results
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
program exchange_test

  use MPI
  use decomp_2d

  implicit none

  integer, parameter :: nx = 8
  integer, parameter :: ny = 8
  integer, parameter :: nz = 8
  integer :: p_row, p_col

  integer :: ierr
  
  ! Initialise
  p_row = 0; p_col = 0
  
  call MPI_Init(ierr)
  if (ierr /= 0) then
     print *, "ERROR: Initialising MPI failed"
  end if
  call decomp_2d_init(nx, ny, nz, p_row, p_col)
  if (nrank == 0) then
     print *, "Testing halo exchange routines"
  end if
  
  ! Run tests
  call test_exchange("X")

  ! Finalise
  call decomp_2d_finalize
  call MPI_Finalize(ierr)
  if (ierr /= 0) then
     print *, "ERROR: Finalising MPI failed"
  end if
  
contains

  !! TEST test_exchange
  !    Tests that the halo exchange in given pencil works.
  !    Given a pre-allocated array with halo buffers this should fill the halo buffers with the
  !    correct data.
  !    Should test a range of halo levels.
  !    - Corner case: zero halo levels.
  subroutine test_exchange(orientation)

    character(len=*), intent(in) :: orientation
    
    real(mytype), dimension(:,:,:), allocatable :: u
    integer :: nlevels
    integer, dimension(3) :: levels

    integer :: xhalo, yhalo, zhalo
    integer, dimension(3) :: starts, ends, sizes

    xhalo = 1; yhalo = 1; zhalo = 1
    if (orientation == "X") then
       xhalo = 0
       starts = xstart
       ends = xend
       sizes = xsize
    else if (orientation == "Y") then
       yhalo = 0
       starts = ystart
       ends = yend
       sizes = ysize
    else if (orientation == "Z") then
       zhalo = 0
       starts = zstart
       ends = zend
       sizes = zsize
    else
       print *, "ERROR: Unknow orientation "//orientation//" test is broken!"
       stop
    end if
    
    do nlevels = 1, 3
       levels=(/ xhalo * nlevels, yhalo * nlevels, zhalo * nlevels /)
       
       allocate(u(1-levels(1):sizes(1)+levels(1), 1-levels(2):sizes(2)+levels(2), 1-levels(3):sizes(3)+levels(3)))

       call local_init(starts, sizes, levels, u)
       if (orientation == "X") then
          call exchange_halo_x(u, opt_xlevel=levels)
       else if (orientation == "Y") then
          call exchange_halo_y(u, opt_ylevel=levels)
       else 
          call exchange_halo_z(u, opt_zlevel=levels)
       end if
       call check(starts, ends, sizes, levels, u)
       
       deallocate(u)
       
    end do
    
  end subroutine test_exchange

  subroutine check(starts, ends, sizes, levels, u)

    integer, dimension(3), intent(in) :: starts, ends, sizes, levels
    real(mytype), dimension(:,:,:), intent(in) :: u
    
    integer :: is, ie, js, je, ks, ke
    integer :: i, j, k
    integer :: io, jo, ko
    integer :: idx

    call valid_range(starts(1), ends(1), sizes(1), nx, levels(1), is, ie)
    call valid_range(starts(2), ends(2), sizes(2), ny, levels(2), js, je)
    call valid_range(starts(3), ends(3), sizes(3), nz, levels(3), ks, ke)

    do k = ks, ke
       ko = k + levels(3)
       do j = js, je
          jo = j + levels(2)
          do i = is, ie
             io = i + levels(1)
             call global_index(i, j, k, starts, idx)
             if (u(io, jo, ko) /= real(idx, mytype)) then
                print *, "ERROR:", nrank, u(io, jo, ko), real(idx, mytype), idx, i, j, k
             end if
          end do
       end do
    end do
  end subroutine check
  
  subroutine local_init(starts, sizes, levels, u)

    integer, dimension(3), intent(in) :: starts, sizes, levels
    real(mytype), dimension(:,:,:), intent(inout) :: u

    integer :: i, j, k, idx
    integer :: io, jo, ko

    u(:,:,:) = 0.0_mytype
    do k = 1, sizes(3)
       ko = k + levels(3)
       do j = 1, sizes(2)
          jo = j + levels(2)
          do i = 1, sizes(1)
             io = i + levels(1)
             call global_index(i, j, k, starts, idx)
             u(io, jo, ko) = real(idx, mytype)
          end do
       end do
    end do
    
  end subroutine local_init
  
  subroutine valid_range(pencil_start, pencil_end, pencil_size, nglobal, nlevels, s, e)

    integer, intent(in) :: pencil_start, pencil_end, pencil_size
    integer, intent(in) :: nglobal
    integer, intent(in) :: nlevels
    integer, intent(out) :: s, e

    if (pencil_start == 1) then
       s = 1
    else
       s = 1 - nlevels
    end if

    if (pencil_end == nglobal) then
       e = pencil_size
    else
       e = pencil_size + nlevels
    end if
    
  end subroutine valid_range
  
  subroutine global_index(i, j, k, starts, idx)

    integer, intent(in) :: i, j, k
    integer, dimension(3), intent(in) :: starts
    integer, intent(out) :: idx

    integer :: ig, jg, kg

    ig = i + starts(1) - 1
    jg = j + starts(2) - 1
    kg = k + starts(3) - 1

    idx = (ig - 1) + (jg - 1) * nx + (kg - 1) * nx * ny
    
  end subroutine global_index
  
end program exchange_test
