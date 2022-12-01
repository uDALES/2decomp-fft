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
  call test_exchange_x

  ! Finalise
  call decomp_2d_finalize
  call MPI_Finalize(ierr)
  if (ierr /= 0) then
     print *, "ERROR: Finalising MPI failed"
  end if
  
contains

  !! TEST test_exchange_x
  !    Tests that the halo exchange in x pencil works.
  !    Given a pre-allocated array with halo buffers this should fill the halo buffers with the
  !    correct data.
  !    Should test a range of halo levels.
  !    - Corner case: zero halo levels.
  subroutine test_exchange_x

    real(mytype), dimension(:,:,:), allocatable :: u
    integer :: nlevels

    integer :: i, j, k
    integer :: is, ie, js, je, ks, ke
    integer :: idx
    
    do nlevels = 1, 3

       allocate(u(xsize(1), 1-nlevels:xsize(2)+nlevels, 1-nlevels:xsize(3)+nlevels))

       u(:,:,:) = 0.0_mytype
       do k = 1, xsize(3)
          do j = 1, xsize(2)
             do i = 1, xsize(1)
                call global_index(i, j, k, xstart, idx)
                u(i, j, k) = real(idx, mytype)
             end do
          end do
       end do

       call exchange_halo_x(u, opt_xlevel=(/ 0, nlevels, nlevels /))

       if (xstart(1) == 1) then
          is = 1
       else
          is = 1 - nlevels
       end if
       if (xend(1) == nx) then
          ie = xsize(1)
       else
          ie = xsize(1) + nlevels
       end if
       if (xstart(2) == 1) then
          js = 1
       else
          js = 1 - nlevels
       end if
       if (xend(2) == ny) then
          je = xsize(2)
       else
          je = xsize(2) + nlevels
       end if
       if (xstart(3) == 1) then
          ks = 1
       else
          ks = 1 - nlevels
       end if
       if (xend(3) == nz) then
          ke = xsize(3)
       else
          ke = xsize(3) + nlevels
       end if
       
       do k = ks, ke
          do j = js, je
             do i = is, ie
                call global_index(i, j, k, xstart, idx)
                if (u(i, j, k) /= real(idx, mytype)) then
                   print *, "ERROR:", nrank, u(i, j, k), real(idx, mytype), idx, i, j, k
                end if
             end do
          end do
       end do
       deallocate(u)
       
    end do
    
  end subroutine test_exchange_x

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
