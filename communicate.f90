!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
!BOP

 module communicate

! !MODULE: communicate
! !DESCRIPTION:
! This module contains the necessary routines and variables for
! communicating between processors.
!
! !REVISION HISTORY:
! CVS:$Id: communicate.F90,v 1.7 2002/05/07 17:39:40 pwjones Exp $
! CVS:$Name: POP_2_0_1 $
!
! !USES:

   use kinds_mod

   implicit none
   private
   save

! !PUBLIC MEMBER FUNCTIONS:
    

   public ::  init_communicate, &
              exit_message_environment, &
              create_comm_group, &
              destroy_comm_group, &
              abort_message_environment, &
	      master_print_message

    interface master_print_message
	module procedure master_print_str
	module procedure master_print_str_dbl
	module procedure master_print_str_int
    end interface 
! !PUBLIC DATA MEMBERS:

   integer (int_kind), public :: &
      mpi_dbl, &! MPI type for dbl_kind
      my_task, group_id, &! MPI task number for this task
      master_task, &! task number of master task
      nprocs

   integer (int_kind), public :: comm, comm_group  ! communicator
   integer (int_kind), parameter, public :: &
      mpitag_bndy_2d = 1, &! MPI tags for various
      mpitag_gs = 1000 ! communication patterns

!EOP
!BOC
!EOC
!***********************************************************************

 contains

!***********************************************************************
!BOP
! !IROUTINE: init_communicate
! !INTERFACE:

 subroutine init_communicate

! !DESCRIPTION:
! This routine sets up MPI environment and defines ocean
! communicator.
!
! !REVISION HISTORY:
! same as module

!EOP
!BOC
!-----------------------------------------------------------------------
!
! local variables
!
!-----------------------------------------------------------------------

   include 'mpif.h' ! MPI Fortran include file

   integer (int_kind) :: ierr ! MPI error flag

!-----------------------------------------------------------------------
!
! initiate mpi environment and create communicator for internal
! ocean communications
!
!-----------------------------------------------------------------------

   call MPI_INIT(ierr)

   master_task = 0
   call MPI_COMM_RANK (MPI_COMM_WORLD, my_task, ierr)
   call MPI_COMM_SIZE (MPI_COMM_WORLD, nprocs, ierr)

!-----------------------------------------------------------------------
!
! On some 64-bit machines where real_kind and dbl_kind are
! identical, the MPI implementation uses MPI_REAL for both.
! In these cases, set MPI_DBL to MPI_REAL.
!
!-----------------------------------------------------------------------
   comm = MPI_COMM_WORLD
   MPI_DBL = MPI_DOUBLE_PRECISION
   !write(*,'(a26, i6, i6)') 'myproc/NPRORC = ', my_task, nprocs

!-----------------------------------------------------------------------
!EOC

 end subroutine init_communicate

!***********************************************************************
!BOP
! !IROUTINE: exit_message_environment
! !INTERFACE:

 subroutine exit_message_environment(ierr)

! !DESCRIPTION:
! This routine exits the message environment properly when model
! stops.
!
! !REVISION HISTORY:
! same as module

! !INCLUDES:

   include 'mpif.h' ! MPI Fortran include file

! !OUTPUT PARAMETERS:

   integer (int_kind), intent(out) :: ierr ! MPI error flag

!EOP
!BOC
!-----------------------------------------------------------------------

   call MPI_FINALIZE(ierr)

!-----------------------------------------------------------------------
!EOC

 end subroutine exit_message_environment

 subroutine create_comm_group(nproc_x, nproc_y, ierr)
 
! !INCLUDES:

   include 'mpif.h' ! MPI Fortran include file
      
! !OUTPUT PARAMETERS:
   integer(int_kind):: &
     MPI_GROUP_OLD, &! group of processors assigned to ocn
     MPI_GROUP_NEW ! group of processors assigned to new dist
  
   integer (int_kind), intent(in)  :: nproc_x, nproc_y
   integer (int_kind), intent(out) :: ierr ! MPI error flag
   integer (int_kind), dimension(3):: range

   integer (int_kind) :: i, rank

   
   !call MPI_COMM_GROUP (MPI_COMM_WORLD, MPI_GROUP_OLD, ierr)
   !do i = 0, nproc_y-1
   !    range(1) = nproc_x * i
   !    range(2) = nproc_x * i + nproc_x-1
   !    range(3) = 1

   !    if (my_task >=range(1) .and. my_task <=range(2) ) then
   !        call MPI_GROUP_RANGE_INCL(MPI_GROUP_OLD, 1, range, &
   !     	MPI_GROUP_NEW, ierr)
   !        write(*,*) "group" , i, my_task, MPI_GROUP_NEW
   !    endif 
   !enddo 

   !call MPI_COMM_CREATE (MPI_COMM_WORLD, MPI_GROUP_NEW, &
   !	comm_group, ierr)

   call MPI_COMM_SPLIT(comm, my_task/nproc_x, mod(my_task, nproc_x),comm_group, ierr)
   !call MPI_COMM_RANK ( comm_group, rank, ierr)
   !write(*,*) "group" ,  my_task, rank, comm_group
   !call MPI_BARRIER(MPI_COMM_WORLD, ierr)


 end subroutine create_comm_group
 subroutine destroy_comm_group()

! !INCLUDES:

   include 'mpif.h' ! MPI Fortran include file

! !OUTPUT PARAMETERS:

   integer (int_kind) :: i,  ierr ! MPI error flag

   !if (allocated(comm_group) ) then 
   !   ! do i = 0, nproc_y-1
   !   !     call MPI_GROUP_FREE ( MPI_GROUP_OLD, ierr)
   !   ! enddo

   !    deallocate(comm_group)
   !endif 



 end subroutine destroy_comm_group

!***********************************************************************
!BOP
! !IROUTINE: abort_message_environment
! !INTERFACE:

 subroutine abort_message_environment(ierr)

! !DESCRIPTION:
! This routine aborts the message environment when model stops.
! It will attempt to abort the entire MPI COMM WORLD.
!
! !REVISION HISTORY:
! same as module

! !INCLUDES:

   include 'mpif.h' ! MPI Fortran include file

! !OUTPUT PARAMETERS:

   integer (int_kind), intent(out) :: ierr ! MPI error flag

!EOP
!BOC
!-----------------------------------------------------------------------

   call MPI_BARRIER(MPI_COMM_WORLD, ierr)
   call MPI_ABORT(MPI_COMM_WORLD, ierr)
   call MPI_FINALIZE(ierr)

!-----------------------------------------------------------------------
!EOC

 end subroutine abort_message_environment
 
 subroutine master_print_str(message)

   include 'mpif.h' ! MPI Fortran include file

   character(*), intent(in) :: message 
   integer (int_kind) :: ierr ! MPI error flag

   if(my_task == master_task) then 
       write(*,*) "===========  ", message, "  =============="
   end if 
       call MPI_BARRIER(MPI_COMM_WORLD, ierr)

 end subroutine master_print_str

 subroutine master_print_str_dbl(num_p, message)

   include 'mpif.h' ! MPI Fortran include file
    
   real(r8), intent(in) :: num_p
   character(*), intent(in) :: message 
   integer (int_kind) :: ierr ! MPI error flag

   call MPI_BARRIER(MPI_COMM_WORLD, ierr)
   if(my_task == master_task) then 
       write(*,*) "_________________________________"
       print *,  message, num_p
       write(*,*) "---------------------------------"
   end if 

 end subroutine master_print_str_dbl

 subroutine master_print_str_int(num_p, message)

   include 'mpif.h' ! MPI Fortran include file
    
   integer(int_kind), intent(in) :: num_p
   character(*), intent(in) :: message 
   integer (int_kind) :: ierr ! MPI error flag

   call MPI_BARRIER(MPI_COMM_WORLD, ierr)
   if(my_task == master_task) then 
       write(*,*) "=========== ", message, num_p, " ==========="
   end if 

 end subroutine master_print_str_int
 
 end module communicate

!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
