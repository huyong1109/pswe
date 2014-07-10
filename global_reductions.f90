!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||

 module global_reductions

!BOP
! !MODULE: global_reductions
! !DESCRIPTION:
!  This module contains all the routines for performing global
!  reductions like global sums, minvals, maxvals, etc.
!
! !REVISION HISTORY:
!  CVS:$Id: global_reductions.F90,v 1.5 2003/01/10 00:11:10 jfd Exp $
!  CVS:$Name: POP_2_0_1 $

! !USES:

   use kinds_mod
   use communicate
   use module_para
   use distribution

   implicit none
   private
   save

! !PUBLIC MEMBER FUNCTIONS:

   public :: global_sum,      &
	     polar_sum,	      &
	     lat_gather,      &
	     lat_scatter,      &
             init_global_reductions

   interface global_sum 
       module procedure global_sum_dbl
       module procedure global_sum_scalar_dbl
   end interface 

 contains

!***********************************************************************
!BOP
! !IROUTINE: init_global_reductions
! !INTERFACE:

 subroutine init_global_reductions

! !DESCRIPTION:
!  Initializes necessary buffers for global reductions.
!
! !REVISION HISTORY:
!  same as module
!EOP
!BOC
!-----------------------------------------------------------------------

  !call get_timer(timer_local, 'SUM_LOCAL')
  !call get_timer(timer_mpi  , 'SUM_MPI')

!-----------------------------------------------------------------------
!EOC

 end subroutine init_global_reductions

!***********************************************************************

 function global_sum_dbl(X)


   include 'mpif.h'  ! MPI Fortran include file

! !INPUT PARAMETERS:

   real (r8), dimension(0:nloc_x+1,0:nloc_y+1), intent(in) :: &
      X                    ! array to be summed

! !OUTPUT PARAMETERS:

   real (r8) :: &
      global_sum_dbl       ! resulting global sum

   real (r8) ::          &
      local_sum           ! sum of all local blocks

   integer (int_kind) :: &
      i,j,             &! local counters
      ierr                ! MPI error flag



   local_sum = .0

   do j =1, nloc_y
       do i =1, nloc_x 
           local_sum = local_sum + X(i,j)
       end do
   end do
   if (jproc == 0 ) then 
       do i=1,nloc_x
          local_sum = local_sum + X(i,0)
       end do
   end if 
   if (jproc == nproc_y -1 ) then 
       do i=1,nloc_x
	  local_sum = local_sum + X(i,nloc_y+1)
       end do
   end if 
   call MPI_ALLREDUCE(local_sum, global_sum_dbl, 1, &
                            mpi_dbl, MPI_SUM, comm, ierr)

 end function global_sum_dbl

!***********************************************************************

 subroutine polar_sum(X, polar_sum_s, polar_sum_n)


   include 'mpif.h'  ! MPI Fortran include file

! !INPUT PARAMETERS:

   real (r8), dimension(0:nloc_x+1,0:nloc_y+1), intent(in) :: &
      X                    ! array to be summed

! !OUTPUT PARAMETERS:

   real (r8) :: &
      polar_sum_s, &       ! resulting global sum
      polar_sum_n

   real (r8) ::          &
      local_sum           ! sum of all local blocks

   integer (int_kind) :: &
      i,j,             &! local counters
      ierr                ! MPI error flag



   local_sum = .0
   
   if (jproc == 0 ) then 
       do i =1, nloc_x 
           local_sum = local_sum + X(i,1)
       end do
       call MPI_ALLREDUCE(local_sum, polar_sum_s, 1, &
                                mpi_dbl, MPI_SUM, comm_group, ierr)
   endif 
   

   local_sum = .0
   if (jproc == nproc_y-1 ) then 
       do i =1, nloc_x 
           local_sum = local_sum + X(i,nloc_y)
       end do
       call MPI_ALLREDUCE(local_sum, polar_sum_n, 1, &
                                mpi_dbl, MPI_SUM, comm_group, ierr)
   endif 
   

 end subroutine polar_sum

 subroutine lat_gather(X, Y)


   include 'mpif.h'  ! MPI Fortran include file

! !INPUT PARAMETERS:

   real (r8), dimension(1:nloc_x*p,1:nloc_y*nproc_x), intent(inout) :: X
   real (r8), dimension(0:nloc_x+1,0:nloc_y+1), intent(in)  :: Y                    ! array to be summed
   real (r8), dimension(1:nloc_x,1:nloc_y)  :: YY                    ! array to be summed


   integer (int_kind) :: &
      i,j,		  &! local counters
      ierr                ! MPI error flag
      YY(:,:) = Y(1:nloc_x, 1:nloc_y)


   call MPI_GATHER(YY(1,1),nloc_x*nloc_y,mpi_dbl, X(1,1),nloc_x*nloc_y, mpi_dbl, 0 ,  comm_group, ierr)  

   end subroutine lat_gather


 subroutine lat_scatter(X, Y)


   include 'mpif.h'  ! MPI Fortran include file

! !INPUT PARAMETERS:

   real (r8), dimension(1:nloc_x,1:nloc_y*nproc_x), intent(in) :: X
   real (r8), dimension(0:nloc_x+1,0:nloc_y+1), intent(inout)  :: Y                    ! array to be summed
   real (r8), dimension(1:nloc_x,1:nloc_y)  :: YY                    ! array to be summed

! !OUTPUT PARAMETERS:


   integer (int_kind) :: &
      i,j,             &! local counters
      ierr                ! MPI error flag

      !write(*,*) 'scatter', my_task,  nloc_x*nloc_y, jproc*nproc_x, jproc
      YY(:,:) = 0.0
      call MPI_SCATTER(X(1,1),nloc_x*nloc_y,mpi_dbl,YY(1,1),nloc_x*nloc_y, mpi_dbl, 0, comm_group, ierr)  
      
      Y(1:nloc_x, 1:nloc_y) = YY(:,:) 
    call  MPI_BARRIER(comm_group, ierr)

 end subroutine lat_scatter

 function global_sum_scalar_dbl(local_scalar)

!-----------------------------------------------------------------------
!
!  this function returns the sum of scalar value across processors
!
!-----------------------------------------------------------------------

   include 'mpif.h'  ! MPI Fortran include file

   real (r8), intent(inout) :: &
      local_scalar                ! local scalar to be compared

   real (r8) :: &
      global_sum_scalar_dbl   ! resulting global sum

   integer (int_kind) :: ierr ! MPI error flag

!-----------------------------------------------------------------------

    call MPI_REDUCE(local_scalar, global_sum_scalar_dbl, 1, &
                            mpi_dbl, 0, MPI_SUM, comm, ierr)

!-----------------------------------------------------------------------

 end function global_sum_scalar_dbl


 end module global_reductions

!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
