!
SUBROUTINE cs
!
   use module_para
   use module_array
   use distribution
!
   implicit none
!
   real*8  :: ai, aj      ! working variables
   integer :: j, ierr		  ! working variable
!
   pi=datan(1d0)*4  
   deta1=4*pi/p*a
   deta=pi/q
   deta2=deta*a*2
   detb=deta*nq
!

   if (my_task == master_task) then
	write(stdout,*) '===', detb,deta, nq
   end if 
   do j=0,nloc_y+1
	  ai    = (j+jglobal)*deta-detb
	  c1(j) = dcos(ai)
	  s1(j) = dsin(ai)
   end do
!
   if (jproc == 0) then 
	ai=1.5*deta-detb
   	c1(0)=dcos(ai)/4
	s1(0)=-1.0d0
   end if 
!   
   if (jproc == nproc_y-1) then 
	ai=(n-0.5)*deta-detb
   	c1(nloc_y+1)=dcos(ai)/4.
	s1(nloc_y+1)=1.0d0
    end if 
! 
!
   do j=0,nloc_y+1
 	  c2(j)=c1(j)
   end do
!
   do j=0,nloc_y+1
  	  ai=s1(j)
	  aj=c1(j)*a
	  f1(j)=2.0d0*omg0*ai
	  f2(j)=ai/aj
   end do
!
   !if (my_task == master_task) then
!	write(stdout,*) c1(:)
   !end if 

   do j=0,nloc_y+1
    
	     !write(*,*)	my_task, "===========", j, "=======", deta1, deta2
 	  ai=deta1*c1(j)
	  aj=deta2*c1(j)
	  c11(j)=1./ai
	  c12(j)=1./aj
	  c13(j)=c11(j)*0.5
	  c14(j)=c12(j)*0.5
	     !write(*,*)	my_task, "===========", j, "======="
   end do
    do j = 0, nprocs
	if (my_task == j ) then 
	     write(*,*)	"===========", j, "======="
   	     write(stdout,'(6f7.3)') c1(:)
	 end if 
	     call  MPI_BARRIER(comm, ierr)
    end do 

!
   return
END SUBROUTINE cs
