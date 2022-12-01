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

  integer, parameter :: nx = 64
  integer, parameter :: ny = 64
  integer, parameter :: nz = 64
  integer :: p_row, p_col

  integer :: ierr
  
  ! Initialise
  p_row = 0; p_col = 0
  
  call MPI_Init(ierr)
  if (ierr /= 0) then
     print *, "ERROR: Initialising MPI failed"
  end if
  call decomp_2d_init(nx, ny, nz, p_row, p_col)
  
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

    print *, "Testing x halo exchange"
    
    do nlevels = 1, 5

       allocate(u(nx, 1-nlevels:ny+nlevels, 1-nlevels:nz+nlevels))
       u(:,:,:) = 0.0_mytype

       call exchange_halo_x(u, opt_xlevel=(/ 0, nlevels, nlevels /))
       
       deallocate(u)
       
    end do
    
  end subroutine test_exchange_x
  
end program exchange_test
