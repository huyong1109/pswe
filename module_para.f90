!
MODULE module_para
!
    use kinds_mod
    use communicate
    use exit_mod
    implicit none
!------------------------constants-------------------- 
    real*8,  parameter       ::    omg0 = 7.292d-5		   ! angular velocity of the earth rotation
    real*8,  parameter       ::    a    = 6371000d0		   ! radius of the earth
!
!------------------------parameters-------------------- 
    real*8                   ::    pi					   !
    real*8                   ::    deta,deta1,deta2,detb   !
!

!------------------------domain ------------------------
    integer		     :: nproc_x, nproc_y


!------------------------grid domain -------------------- 
    integer		     ::    p  			   ! to define zonal resolution
    integer		     ::    kn 			   ! to define zonal resolution
    integer		     ::    np 			   ! zonal grid nmuber 
    integer		     ::    n1 		   !
!	   		              							       !
    integer		     ::    q0 		   ! to define meridional resolution
    integer		     ::    nq 		   !
    integer		     ::    q  		   !
    integer		     ::    n  			   ! meridional grid number
! 														   !

    integer	             ::    nm	! Total grid number
!
!------------------------Time management -------------------- 
    integer       ::    tlo              ! time stepsize
    integer       ::    t0               ! initial time
    integer       ::    t1            ! final time
    													   !
    integer       ::    n0 = 0                  ! 0: using RH waves as initial condition (IC); 1: read IC from files
    integer       ::    nyn = 0                 ! Screen output 0: at each step; 1: once 12 hours
    integer       ::    thalf = 432000          ! The time for saving the intermediate result
!
!   When n0 = 1
!														   !
!---------------------- I/O -----------------------------
    integer (i4), parameter, public :: &
	    nml_in = 10, & 
	    stdin = 5,   &
	    stdout = 6

    character*7, parameter, public :: &
   	    nml_filename = 'pswe_in'
   
    character*7, parameter   ::   &
        fu='uui.dat', &		       ! initial field of zonal wind for self-defined experiment
        fv='vvi.dat', &	   	       ! initial field of meridional wind for self-defined experiment
        fh='hhi.dat'	 	       ! initial field of geopotential height for self-defined experiment
   !														   !
    character*6, parameter   ::    &
        ffu='ui.dat', & 		   ! result of zonal wind for self-defined experiment
        ffv='vi.dat', & 		   ! result of meridional wind for self-defined experiment
        ffh='hi.dat' 		   ! result of geopotential height for self-defined experiment


	public :: init_para 
	
contains 

 subroutine init_para

!----------------------------------------------------------------------
!
! local variables
!
!----------------------------------------------------------------------

   integer (int_kind) :: &
      nml_error, sigAbort ! namelist read error flag

!----------------------------------------------------------------------
!
! input namelists
!
!----------------------------------------------------------------------

   namelist /dist_nml/  nproc_x, nproc_y
   namelist /domain_nml/  p, q
   namelist /time_manager_nml/  tlo, t0, t1

!----------------------------------------------------------------------
!
! read domain information from namelist input
!
!----------------------------------------------------------------------

   p = -1
   q = -1
   tlo = -1
   t0 = -1
   t1 = -1

   if (my_task == master_task) then
      open (nml_in, file=nml_filename, status='old',iostat=nml_error)
      if (nml_error /= 0) then
         nml_error = -1
      else
         nml_error = 1
      endif
      do while (nml_error > 0)
         read(nml_in, nml=dist_nml,iostat=nml_error)
      end do
      do while (nml_error > 0)
	  read(nml_in, nml=domain_nml,iostat=nml_error)
      end do
      do while (nml_error > 0)
         read(nml_in, nml=time_manager_nml,iostat=nml_error)
      end do
      if (nml_error == 0) close(nml_in)

      !call broadcast_scalar(nml_error, master_task)
      if (nml_error /= 0) then
         call exit_PSWE(sigAbort,'ERROR reading pswe_in')
      endif

   !call broadcast_scalar(nproc_x, master_task)
   !call broadcast_scalar(nproc_y, master_task)
   !call broadcast_scalar(tlo, master_task)
   !call broadcast_scalar(t0, master_task)
   !call broadcast_scalar(t1, master_task)

!---------------------------------------------------------------------- !
! perform some basic checks on domain
!
!----------------------------------------------------------------------

   !if (my_task == master_task) then
     write(stdout,'(a18)') 'Domain Information'
     write(stdout,'(a26,i6)') '    proc_x = ', nproc_x
     write(stdout,'(a26,i6)') '    proc_y = ', nproc_y
     write(stdout,'(a26,i6)') '    nx_global = ', p
     write(stdout,'(a26,i6)') '    ny_global = ', q
     write(stdout,'(a18)') 'Time   Information'
     write(stdout,'(a26,i6)') '    Start step = ', t0
     write(stdout,'(a26,i6)') '    End   step = ', t1
     write(stdout,'(a26,i6)') '    Step  size = ', tlo

   !endif

   if (p < 1 .or. q  < 1 ) then
      !***
      !*** domain size zero or negative
      !***
      call exit_PSWE(sigAbort,'Invalid domain: size < 1') ! no domain
   endif

    
    kn = p/2				   ! to define zonal resolution
    np = p+1				   ! zonal grid nmuber 
    n1 = np+1			   !

    q0 = q/2
    nq = q0+1			   !
    n  = q+1				   ! meridional grid number
! 														   !

    nm = (np-1)*n		   ! Total grid number
!
     										   !

   !if (my_task == master_task) then
     write(stdout,'(a18)') 'Domain Information'
     write(stdout,'(a26,i6)') '    proc_x = ', nproc_x
     write(stdout,'(a26,i6)') '    proc_y = ', nproc_y
     write(stdout,'(a26,i6)') '    nx_global = ', p
     write(stdout,'(a26,i6)') '    ny_global = ', q
     write(stdout,'(a18)') 'Time   Information'
     write(stdout,'(a26,i6)') '    Start step = ', t0
     write(stdout,'(a26,i6)') '    End   step = ', t1
     write(stdout,'(a26,i6)') '    Step  size = ', tlo

   !endif

   endif

 end subroutine init_para

END MODULE module_para
!
