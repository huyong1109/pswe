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

      integer (int_kind), public :: &
	 iproc,	  &
	 jproc,	  &
	 iglobal, &
	 jglobal, &
	 nloc_x,  &
	 nloc_y
      integer (int_kind), dimension(:,:), allocatable, public :: dist


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

   !include 'mpif.h' ! MPI Fortran include file
   use mpi

   integer (int_kind), intent(in) :: &
	nprocs,		&
   	nproc_x, nproc_y, &  ! number of processors in this distribution
   	nglob_x, nglob_y	! number of processors in this distribution

    ! local parameter 
    integer (int_kind) :: &
	i,j,k,  ierr, 	&
	tmpx,	&
	tmpy
	
    !write(stdout,'(a6,i4, i4, i6, i6, i4, i4, i4 )') '1id',  &
!	& my_task, iproc, jproc, iglobal, jglobal, nloc_x, nloc_y 
    if (.not. allocated( dist) ) allocate(dist(0:nproc_x+1, 0:nproc_y+1))
    dist(:,:) = MPI_PROC_NULL
    call  MPI_BARRIER(comm, ierr)
    
    tmpx  = nglob_x/nproc_x
    tmpy  = (nglob_y-1)/nproc_y

    do i = 0, nproc_x-1
	do j = 0, nproc_y-1
	    k = j*nproc_x+i
	    dist(i+1, j+1) = k
	   
	   if (my_task == k ) then 
	    iproc = i
	    jproc = j
	    iglobal = iproc*tmpx 
	    jglobal = jproc*tmpy +1
	    nloc_x  = tmpx
	    !write(*,*)  my_task,  nglob_x -iglobal, tmpx, nloc_x
	    if (jproc == nproc_y-1) then 
		nloc_y  =nglob_y-jglobal-1
	    else 
		nloc_y = tmpy
	    endif 

	    end if 
	end do 
    end do 

    do j = 1, nproc_y
        dist(0, j) = dist(nproc_x,j)
        dist(nproc_x+1, j) = dist(1,j)
    end do 

    !nloc_x = nloc_x !+ 1  ! grid range x ~ [1, nloc_x]
    !nloc_y = nloc_y !+ 1	 ! grid range y ~ [1, nloc_y]
    !if (my_task == master_task) then
    !    !@write(*,*)	"======================================================="
    !endif 
    !do k = 0, nprocs
    !    if (my_task == k ) then 
    !     !    write(stdout,'(a4,i4, i4, i4, i6, i6, i4, i4 )') 'id2',  &
    !     !        & my_task, iproc, jproc, iglobal, jglobal, nloc_x, nloc_y 
    !         !do i = 0, nproc_x+1
    !         !   write(*,*) dist(i,:)
    !         !end do 
    !     end if 
    !         call  MPI_BARRIER(comm, ierr)
    !end do 

    end subroutine  create_distribution
   

    ! warning!!!!!  please call this sub somewhere 
    subroutine  destroy_distribution 

	if (allocated (dist)) deallocate (dist)
    
    end subroutine  destroy_distribution





end module distribution

!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
