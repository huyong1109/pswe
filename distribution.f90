!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||

 module distribution

!BOP
! !MODULE: distribution
!
! !DESCRIPTION:
! This module provides data types and routines for distributing
! blocks across processors.
!
! !REVISION HISTORY:
! CVS:$Id: distribution.F90,v 1.11 2003/12/23 22:11:40 pwjones Exp $
! CVS:$Name: POP_2_0_1 $

! !USES:

   use module_para
   use communicate

   implicit none
   private
   save

! !PUBLIC TYPES:

      integer (int_kind) :: &
	 iproc,	  &
	 jproc,	  &
	 iglobal, &
	 jglobal, &
	 nloc_x,  &
	 nloc_y


! !PUBLIC MEMBER FUNCTIONS:

   public :: create_distribution

!EOP
!BOC
!EOC
!***********************************************************************

 contains

!***********************************************************************
!BOP
! !IROUTINE: create_distribution
! !INTERFACE:

 subroutine create_distribution(nprocs, nproc_x, nproc_y, nglob_x, nglob_y)

! !INPUT PARAMETERS:


   integer (int_kind), intent(in) :: &
	nprocs,		&
   	nproc_x, nproc_y, &  ! number of processors in this distribution
   	nglob_x, nglob_y	! number of processors in this distribution

    ! local parameter 
    integer (int_kind) :: &
	i,	&
	tmpx,	&
	tmpy	&
	
      
    do i = 0, nprocs
	if my_task == i 
	    iproc = i%nproc_x   
	    jproc = i/nproc_x
	    tmpx  = nglob_x/iproc
	    tmpy  = nglob_y/jproc
	    iglobal = iproc*tmpx
	    jglobal = jproc*tmpy
	    nloc_x  = max(nglob_x-iglobal,tmpx) 
	    nloc_y  = max(nglob_y-jglobal,tmpy) 
	end if 
    end do 
    


 subroutine  create_distribution




end module distribution

!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
