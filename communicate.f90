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

   public :: init_communicate, &
              exit_message_environment, &
              abort_message_environment, &
              get_num_procs, &

! !PUBLIC DATA MEMBERS:

   integer (int_kind), public :: &
      mpi_dbl, &! MPI type for dbl_kind
      my_task, &! MPI task number for this task
      master_task, &! task number of master task
      nprocs

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

!-----------------------------------------------------------------------
!
! On some 64-bit machines where real_kind and dbl_kind are
! identical, the MPI implementation uses MPI_REAL for both.
! In these cases, set MPI_DBL to MPI_REAL.
!
!-----------------------------------------------------------------------

   MPI_DBL = MPI_DOUBLE_PRECISION

!-----------------------------------------------------------------------
!EOC

 end subroutine init_communicate

!***********************************************************************
!BOP
! !IROUTINE: get_num_procs
! !INTERFACE:

 function get_num_procs()

! !DESCRIPTION:
! This function returns the number of processor assigned to
! MPI_COMM_OCN
!
! !REVISION HISTORY:
! same as module

! !OUTPUT PARAMETERS:

   integer (int_kind) :: get_num_procs

!EOP
!BOC
!-----------------------------------------------------------------------
!
! local variables
!
!-----------------------------------------------------------------------

   integer (int_kind) :: ierr

!-----------------------------------------------------------------------

   call MPI_COMM_SIZE(MPI_COMM_OCN, get_num_procs, ierr)

!-----------------------------------------------------------------------
!EOC

 end function get_num_procs

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

   call MPI_BARRIER(MPI_COMM_OCN, ierr)
   call MPI_ABORT(MPI_COMM_WORLD, ierr)
   call MPI_FINALIZE(ierr)

!-----------------------------------------------------------------------
!EOC

 end subroutine abort_message_environment


 end module communicate

!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
