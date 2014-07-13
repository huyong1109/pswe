!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||

 module boundary

! !DESCRIPTION:
!  This module contains data types and routines for updating ghost cell
!  boundaries using MPI calls
!
! !REVISION HISTORY:
!  CVS:$Id: boundary.F90,v 1.13 2004/01/07 19:56:32 pwjones Exp $
!  CVS:$Name: POP_2_0_1 $

! !USES:

   use kinds_mod
   use communicate
   use module_para
   use distribution
   use exit_mod
   !use timers

   implicit none
   private
   save

! !PUBLIC TYPES:

   integer (int_kind), public :: &
       e_proc        ,&! dest   proc for east-west send message
       w_proc        ,&! source proc for east-west recv message
       n_proc        ,&! dest   proc for north-south send message
       s_proc          ! source proc for north-south recv message


! !PUBLIC MEMBER FUNCTIONS:

   public :: create_boundary,  &
             destroy_boundary, &
             update_boundary, &
	     update_latitude, &
	     update_neigbor 

   interface update_boundary  ! generic interface
      module procedure boundary_2d_dbl
   end interface

!EOP
!BOC

contains

!***********************************************************************
!BOP
! !IROUTINE: create_boundary
! !INTERFACE:

 subroutine create_boundary
    
   integer :: i,j, ierr

   do j=1,nproc_y
	do i=1,nproc_x 
	if (dist(i,j)  == my_task) then
	    e_proc = dist(i+1, j)
	    w_proc = dist(i-1, j)
	    n_proc = dist(i, j+1)
	    s_proc = dist(i, j-1)
  
	endif
	end do 
   end do


 end subroutine create_boundary

!***********************************************************************
!BOP
! !IROUTINE: destroy_boundary
! !INTERFACE:

 subroutine destroy_boundary

! !INPUT/OUTPUT PARAMETERS:


 end subroutine destroy_boundary

!***********************************************************************
!BOP
! !IROUTINE: update_boundary
! !INTERFACE:

 subroutine boundary_2d_dbl(ARRAY)

   include 'mpif.h'   ! MPI Fortran include file

! !INPUT PARAMETERS:
   real (r8), dimension(0:nloc_x+1,0:nloc_y+1), intent(inout) :: &
      ARRAY              ! array containing horizontal slab to update

!-----------------------------------------------------------------------
!
!  local variables
!
!-----------------------------------------------------------------------

   integer (int_kind) ::           &
      i,j,    	                   &! dummy loop indices
      ierr                          ! MPI error flag

   integer (int_kind) ::	&
      snd_request,              &! MPI request ids
      rcv_request                ! MPI request ids

   integer (int_kind) :: &
      snd_status,               &! MPI status flags
      rcv_status                 ! MPI status flags

   real (r8), dimension(:), allocatable :: &
      buf_ew_snd,       &! message buffer for east-west sends
      buf_ew_rcv,       &! message buffer for east-west recvs
      buf_ns_snd,       &! message buffer for north-south sends
      buf_ns_rcv         ! message buffer for north-south recvs


!-----------------------------------------------------------------------
!
!  allocate buffers for east-west sends and receives
!
!-----------------------------------------------------------------------
   !call timer_start(bndy_2d_recv)
   allocate (buf_ew_snd(1:nloc_y), &
	     buf_ew_rcv(1:nloc_y), &
	     buf_ns_snd(1:nloc_x), &
	     buf_ns_rcv(1:nloc_x))
	     
   !==================send to east, rcv from west ============

   call MPI_IRECV(buf_ew_rcv(1), nloc_y, mpi_dbl,   &
                     w_proc, 1, comm, rcv_request, ierr)
   
   if (e_proc /= MPI_PROC_NULL) then 
       do  j = 1, nloc_y
	   buf_ew_snd(j)  = ARRAY(nloc_x, j)
       end do 
   end if 

   call MPI_ISEND(buf_ew_snd(1), nloc_y, mpi_dbl, &
                     e_proc, 1, comm, snd_request, ierr)

   call MPI_WAIT(rcv_request, rcv_status, ierr)
   

   if (w_proc /= MPI_PROC_NULL) then 
       do  j = 1, nloc_y
	   ARRAY(0, j) = buf_ew_rcv(j)
       end do 
   end if

   call MPI_WAIT(snd_request, snd_status, ierr)

   !==================send to west, rcv from east============
   buf_ew_rcv(:) = 0.
   call MPI_IRECV(buf_ew_rcv(1), nloc_y, mpi_dbl,   &
                     e_proc, 2, comm, rcv_request, ierr)

   if (w_proc /= MPI_PROC_NULL) then 
   do  j = 1, nloc_y
       buf_ew_snd(j)  = ARRAY(1, j)
   end do 
   end if

   call MPI_ISEND(buf_ew_snd(1), nloc_y, mpi_dbl, &
                     w_proc, 2, comm, snd_request, ierr)

   call MPI_WAIT(rcv_request, rcv_status, ierr)

   if (e_proc /= MPI_PROC_NULL) then 
   do  j = 1, nloc_y
       ARRAY(nloc_x+1, j) = buf_ew_rcv(j)
   end do 
   end if 

   call MPI_WAIT(snd_request, snd_status, ierr)
    

   deallocate(buf_ew_snd, buf_ew_rcv)

   !==================send to north, rcv from south ============

   call MPI_IRECV(buf_ns_rcv(1), nloc_x, mpi_dbl,   &
                     s_proc, 3, comm, rcv_request, ierr)

   if (n_proc /= MPI_PROC_NULL) then 
   do  i = 1, nloc_x
       buf_ns_snd(i)  = ARRAY(i, nloc_y)
   end do
   end if 

   call MPI_ISEND(buf_ns_snd(1), nloc_x, mpi_dbl, &
                     n_proc, 3, comm, snd_request, ierr)

   call MPI_WAIT(rcv_request, rcv_status, ierr)

   if (s_proc /= MPI_PROC_NULL) then 
   do  i = 1, nloc_x
       ARRAY(i, 0) = buf_ns_rcv(i)
   end do 
   end if 

   call MPI_WAIT(snd_request, snd_status, ierr)

   !==================send to south, rcv from north ============
   buf_ns_rcv(:) = 0.
   call MPI_IRECV(buf_ns_rcv(1), nloc_x, mpi_dbl,   &
                     n_proc, 4, comm, rcv_request, ierr)

   if (s_proc /= MPI_PROC_NULL) then 
   do  i = 1, nloc_x
       buf_ns_snd(i)  = ARRAY(i, 1)
   end do 
   end if 

   call MPI_ISEND(buf_ns_snd(1), nloc_x, mpi_dbl, &
                     s_proc, 4, comm, snd_request, ierr)

   call MPI_WAIT(rcv_request, rcv_status, ierr)

   if (n_proc /= MPI_PROC_NULL) then 
   do  i = 1, nloc_x
       ARRAY(i, nloc_y+1) = buf_ns_rcv(i)
   end do 
   end if 

   call MPI_WAIT(snd_request, snd_status, ierr)

   deallocate(buf_ns_snd, buf_ns_rcv)


 end subroutine boundary_2d_dbl

 subroutine  update_latitude(ARRAY)
   include 'mpif.h'   ! MPI Fortran include file

! !INPUT PARAMETERS:
   real (r8), dimension(0:nloc_x+1,0:nloc_y+1), intent(inout) :: &
      ARRAY              ! array containing horizontal slab to update

!-----------------------------------------------------------------------
!
!  local variables
!
!-----------------------------------------------------------------------

   integer (int_kind) ::           &
      i,j,    	                   &! dummy loop indices
      ierr                          ! MPI error flag

   integer (int_kind) ::	&
      snd_request,              &! MPI request ids
      rcv_request                ! MPI request ids

   integer (int_kind) :: &
      snd_status,               &! MPI status flags
      rcv_status                 ! MPI status flags

   real (r8), dimension(:), allocatable :: &
      buf_ew_snd,       &! message buffer for east-west sends
      buf_ew_rcv,       &! message buffer for east-west recvs
      buf_ns_snd,       &! message buffer for north-south sends
      buf_ns_rcv         ! message buffer for north-south recvs


!-----------------------------------------------------------------------
!
!  allocate buffers for east-west sends and receives
!
!-----------------------------------------------------------------------
   !call timer_start(bndy_2d_recv)
   allocate (buf_ew_snd(1:nloc_y), &
	     buf_ew_rcv(1:nloc_y), &
	     buf_ns_snd(1:nloc_x), &
	     buf_ns_rcv(1:nloc_x))
	     
   !==================send to east, rcv from west ============

   call MPI_IRECV(buf_ew_rcv(1), nloc_y, mpi_dbl,   &
                     w_proc, 1, comm, rcv_request, ierr)
   
   if (e_proc /= MPI_PROC_NULL) then 
       do  j = 1, nloc_y
	   buf_ew_snd(j)  = ARRAY(nloc_x, j)
       end do 
   end if 

   call MPI_ISEND(buf_ew_snd(1), nloc_y, mpi_dbl, &
                     e_proc, 1, comm, snd_request, ierr)

   call MPI_WAIT(rcv_request, rcv_status, ierr)
   

   if (w_proc /= MPI_PROC_NULL) then 
       do  j = 1, nloc_y
	   ARRAY(0, j) = buf_ew_rcv(j)
       end do 
   end if

   call MPI_WAIT(snd_request, snd_status, ierr)

   !==================send to west, rcv from east============
   buf_ew_rcv(:) = 0.
   call MPI_IRECV(buf_ew_rcv(1), nloc_y, mpi_dbl,   &
                     e_proc, 2, comm, rcv_request, ierr)

   if (w_proc /= MPI_PROC_NULL) then 
   do  j = 1, nloc_y
       buf_ew_snd(j)  = ARRAY(1, j)
   end do 
   end if

   call MPI_ISEND(buf_ew_snd(1), nloc_y, mpi_dbl, &
                     w_proc, 2, comm, snd_request, ierr)

   call MPI_WAIT(rcv_request, rcv_status, ierr)

   if (e_proc /= MPI_PROC_NULL) then 
   do  j = 1, nloc_y
       ARRAY(nloc_x+1, j) = buf_ew_rcv(j)
   end do 
   end if 

   call MPI_WAIT(snd_request, snd_status, ierr)
    

   deallocate(buf_ew_snd, buf_ew_rcv)


 end subroutine update_latitude
 
 subroutine  update_neigbor 
 

 end subroutine update_neigbor

end module boundary

!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
