#define DAT_FNAME ("mesh240_orig.txt")

#define CACHE_SIZE_L2 (1 * 1024 * 1024)

Module vars
Implicit None
save

  ! *** parameters to control iteration counts ***
  integer, parameter :: nouter  = 20, & ! number of times to loop through all kernels (running each ninner times)
                        ninner  = 1, & ! number of times to repeat each kernel successively
                        nkernel = 1  ! number of kernels run below (used in storing times, solutions)
  ! --
  ! profiling variables
  integer :: it1, it2
  integer :: m, mz, n, ndiag, m_zz
  double precision, dimension(nouter) :: wt
  double precision :: wtb, wtk, std, mean, median

  ! problem variables
  integer :: max_nEdgesOnEdge, nEdges, nCells
  integer, allocatable, dimension(:) :: nEdgesOnEdge
  integer, allocatable, dimension(:,:) :: cellsOnEdge, edgesOnEdge
  double precision, allocatable, dimension(:) :: fEdge,normalBarotropicVelocitySubcycleCur
  double precision, allocatable, dimension(:,:) :: weightsOnEdge
  double precision, allocatable, dimension(:) :: sshSubcycleCur
  double precision, allocatable, dimension(:) :: dcEdge
  double precision, allocatable, dimension(:) :: barotropicForcing
  double precision, allocatable, dimension(:) :: normalBarotropicVelocitySubcycleNew
  integer, allocatable, dimension(:) :: edgeMask

  ! --
  integer :: i, j, k, mb, d, num_threads=1, kernel=-1, s
  integer :: numtasks, rank, ierr
  double precision :: norm2

  ! miscellaneous
  double precision, allocatable, dimension(:) :: cache_buffer
  integer :: buffer_size = CACHE_SIZE_L2 * 64

  character(len=32) :: args
  integer :: argcount


contains 

subroutine init_vars()
!  argcount = command_argument_count()
!  if (argcount.ge.2) then
!    call get_command_argument(2, args)
!    read(args,*)nRHS
!    write(*,*)'Using RHS block size: ',nRHS
!  else
!    write(*,*)'Using default RHS block size: ',nRHS
!  endif
end subroutine init_vars

end module vars




program kernels
#if defined(USE_SDE)
  use itt_sde_fortran
#endif
  use vars

  implicit none
  include 'mpif.h'
  include 'omp_lib.h'

  print *, "Cache cleaning buffer size = ", buffer_size * 8 / 1024 / 1024, "MB"
  allocate(cache_buffer(buffer_size))
  call random_seed()
  do i = 1, buffer_size
    call random_number(cache_buffer(i))
  enddo
  
  ! use bulk vector mapping instead of inlined one
  ! --
  call MPI_Init(ierr)
  call MPI_Comm_size(MPI_COMM_WORLD, numtasks, ierr)
  call MPI_Comm_rank(MPI_COMM_WORLD, rank, ierr)
  
  !initialize variables
!  call init_vars()

  ! load the data
  ! here, the data are packed such that each SIMD height x ndiag width panel of
  ! the matrix element and indexing data has been stored contiguously
  print *, "starting to read initial conditions"
  open(unit=10, file=DAT_FNAME, action='read')
  read(10, '(i10)') max_nEdgesOnEdge
  read(10, '(i10)') nEdges
  read(10, '(i10)') nCells

  print *, max_nEdgesOnEdge, nEdges, nCells

  allocate(cellsOnEdge(2,nEdges))
  allocate(nEdgesOnEdge(nEdges))
  allocate(edgesOnEdge(max_nEdgesOnEdge, nEdges))
  allocate(weightsOnEdge(max_nEdgesOnEdge, nEdges))
  allocate(normalBarotropicVelocitySubcycleCur(nEdges))
  allocate(fEdge(nEdges))
  allocate(sshSubcycleCur(nCells))
  allocate(dcEdge(nEdges))
  allocate(barotropicForcing(nEdges))
  allocate(edgeMask(nEdges))
  allocate(normalBarotropicVelocitySubcycleNew(nEdges))


  read(10, '(i10)' ) cellsOnEdge
  read(10, '(i10)' ) edgesOnEdge
  read(10, '(i10)' ) nEdgesOnEdge
  read(10, '(e23.16)') weightsOnEdge
  read(10, '(e23.16)') normalBarotropicVelocitySubcycleCur
  read(10, '(e23.16)') fEdge
  read(10, '(e23.16)') sshSubcycleCur
  read(10, '(e23.16)') dcEdge
  read(10, '(e23.16)') barotropicForcing
  close(10)

  print *, 'Done reading data'

  do i = 1, nEdges
    edgeMask(i) = 1
  enddo


  ! run once before timing and before initialization to get OMP initialized and for first touch policy
  !$omp parallel
  call run_kernel()
  !$omp end parallel

  ! ... run timing tests ...

  !call MPI_Pcontrol(1, "Measure BW"//char(0))

  print *, "starting iterations..."
  do it1 = 1, nouter
    call clear_cache()
    wtb = MPI_Wtime()
    do it2 = 1, ninner
#if defined(USE_SDE) || defined(USE_VTUNE)
      call start_collection()
#endif
      !$omp parallel
      call run_kernel()
      !$omp end parallel
#if defined(USE_SDE) || defined(USE_VTUNE)
      call stop_collection()
#endif
    end do ! it2
    wt(it1) = MPI_Wtime() - wtb
    print *, "Iteration", it1, "time: ", wt(it1)
  end do ! it1

  !call MPI_Pcontrol(-1, "Measure BW"//char(0))

  wtk = sum(wt)
  mean = wtk / nouter / ninner
  call sort(wt, size(wt))
  median = wt(nouter / 2) / ninner
  std = SQRT(sum((wt(1:nouter) / ninner - mean) ** 2) / nouter)
  write(*,*) 'Time: total = ', sum(wt), ' median = ', median
  write(*,*) '        min = ', minval(wt) / ninner, '    max = ', maxval(wt) / ninner
  write(*,*) '       mean = ', mean, ' stddev = ', std

  !$omp parallel
  num_threads = omp_get_num_threads()
  !$omp do schedule(runtime) reduction(+:norm2)
  do i = 1, nEdges
    norm2 = norm2 + normalBarotropicVelocitySubcycleNew(i) * normalBarotropicVelocitySubcycleNew(i)
  enddo
  !$omp end do
  !$omp end parallel
  print *, "Number of OpenMP threads:", num_threads
  print *, "Norm value:", sqrt(norm2)

  call MPI_Finalize(ierr)
end program kernels


subroutine clear_cache()
  use vars
  implicit none
  double precision :: rsum
  do i = 1, buffer_size
    call random_number(cache_buffer(i))
  enddo
  rsum = 0.0
  !$omp do schedule(runtime) reduction(+:rsum)
  do i = 1, buffer_size
    rsum = rsum + cache_buffer(i)
  enddo
  !$omp end do
  rsum = rsum / buffer_size
  print *, "Mean random number: ", rsum
end subroutine clear_cache


subroutine run_kernel()
  use vars
  implicit none

#ifdef __KERNEL_1
  kernel = 1
  call edge_bench(nEdges, nCells, max_nEdgesOnEdge, nEdgesOnEdge, cellsOnEdge, &
                  edgesOnEdge, fEdge, normalBarotropicVelocitySubcycleCur, &
                  weightsOnEdge, sshSubcycleCur, dcEdge, barotropicForcing, edgeMask, &
                  normalBarotropicVelocitySubcycleNew)
#endif

end subroutine run_kernel


! kernel 1
subroutine edge_bench(nEdges, nCells, max_nEdgesOnEdge, nEdgesOnEdge, cellsOnEdge, &
                   edgesOnEdge, fEdge, normalBarotropicVelocitySubcycleCur, &
                   weightsOnEdge, sshSubcycleCur, dcEdge, barotropicForcing, edgeMask, &
                   normalBarotropicVelocitySubcycleNew)
  implicit none

  ! arguments
  integer :: nEdges, nCells, max_nEdgesOnEdge
  integer, dimension(nEdges) :: nEdgesOnEdge
  integer, dimension(2, nEdges) :: cellsOnEdge
  integer, dimension(max_nEdgesOnEdge, nEdges) :: edgesOnEdge
  double precision, dimension(nEdges) :: fEdge
  double precision, dimension(nEdges) :: normalBarotropicVelocitySubcycleCur
  double precision, dimension(max_nEdgesOnEdge, nEdges) :: weightsOnEdge
  double precision, dimension(nCells) :: sshSubcycleCur
  double precision, dimension(nEdges) :: dcEdge
  double precision, dimension(nEdges) :: barotropicForcing
  integer, dimension(nEdges) :: edgeMask
  double precision, dimension(nEdges) :: normalBarotropicVelocitySubcycleNew

  ! local
  double precision :: CoriolisTerm
  integer :: iEdge, i, eoe, cell1, cell2
  double precision :: const1, const2
  const1 = 360.0
  const2 = 9.86

  !$omp do schedule(runtime) private(cell1, cell2, CoriolisTerm, i, eoe)
  do iEdge = 1, nEdges

    cell1 = cellsOnEdge(1,iEdge)
    cell2 = cellsOnEdge(2,iEdge)

    ! Compute the barotropic Coriolis term, -f*uPerp
    CoriolisTerm = 0.0
!    !dir$ VECTOR ALIGNED
    do i = 1, nEdgesOnEdge(iEdge)
       eoe = edgesOnEdge(i,iEdge)
       CoriolisTerm = CoriolisTerm + weightsOnEdge(i,iEdge) &
                      * normalBarotropicVelocitySubcycleCur(eoe) * fEdge(eoe)
    end do

    normalBarotropicVelocitySubcycleNew(iEdge) &
      = (normalBarotropicVelocitySubcycleCur(iEdge) &
        + const1 * (CoriolisTerm - const2 &
        * (sshSubcycleCur(cell2) - sshSubcycleCur(cell1) ) &
        / dcEdge(iEdge) + barotropicForcing(iEdge))) * edgeMask(iEdge)

  end do
  !$omp end do

end subroutine edge_bench

subroutine  sort(x, s)
  implicit none
  integer :: s
  double precision, dimension(s) :: x

  integer :: i, ii, loc
  double precision :: minimum, tmp
  do i = 1, s-1
    minimum  = x(i)
    loc = i
    do ii = i+1, s
      if (x(ii) < minimum) then
        minimum  = x(ii)
        loc = ii
      endif
    enddo
    tmp  = x(i)
    x(i) = x(loc)
    x(loc) = tmp
  enddo
end subroutine sort

