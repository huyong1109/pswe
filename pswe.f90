!*****************************************************************************!
!				Parallel Version of                           !
!			       The Barotropic Model on Arakawa A-grid	      !
!                               by Bin Wang	                              !
!	        with the implicit scheme of energy conservation               !
!*****************************************************************************!	
!
program PSWE
    use kinds_mod
    use communicate
    use module_para
    use module_array
    use distribution
    use global_reductions
    use boundary
    use module_io

    implicit none 
    ! local 
    integer(int_kind) :: i, j 
    integer :: ierr
    real(r8) :: global_s
    
    real(r8) :: tener,tener0	! total energy at tn and t0, respectively
    real(r8) :: tmass,tmass0	! total mass at tn and t0, respectively
!																  !
    real(r8), dimension(:,:), allocatable :: pu,pv,ph                       ! for GrADS saving
!																  !
    real(r8) :: ai,dt                          ! working variables
!																  !
    integer  :: tlp,iter,irecu,irecv,irech     ! working variables
    integer  :: iws,nt,nw,iwr,tt               !	working variables
!																  !
    real(r8), external            :: inner, mass						  ! a external function to calculate inner product
!
    call  init_communicate
    call  master_print_message("Init Parameter")
    call  init_para
    call  master_print_message("Init Array")
    call  init_array
    call  master_print_message("Create Boundary")
    call  create_boundary 
    call  create_comm_group(nproc_x, nproc_y, ierr)
    
    tt = t0  
    iwr=(t1-tt)/tlo	!hy total step 
    tlp=t1-tt-tlo*iwr	!hy last timestep
    
    if (tlp.gt.0) then
       iwr=iwr+1
    else
       tlp=tlo
    end if
!
    if (my_task == master_task ) then
	print *,'Initial time is',tt
	print *,'Final time is',t1
	print *,'Time stepsize is',tlo

	if (n0.eq.0) then
	   print *,'This is a RH-wave experiment,'
	   print *,'i.e., the IC is RH-wave'
	else
	   print *,'This is self-defined experiment,'
	   print *,'i.e., the IC is read from files'
	end if
!
	if (nyn.ne.0) then
	   print *,'The Energy... will be shown once 12 hours'
	else
	   print *,'The Energy... will be shown at each time step'
	end if
    end if
   
    
    
    
    call  master_print_message("Compute parameters")
    call  cs  
    
    !      Using Rossby-Haurwitz waves as initial condition
    call  haurwitz  
    ! variables transform
    wu(:,:) = 0.0
    wv(:,:) = 0.0
    do j=1,nloc_y
	do i=1,nloc_x
	    ai=dsqrt(wh(i,j))
	    wu(i,j)=u(i,j)*ai
	    wv(i,j)=v(i,j)*ai
	end do
    end do
    
    !call mpi_print(wv) 
    
    tener0=inner(wu,wv,wh,wu,wv,wh)
    tmass0 = mass(wh)


    nt=23
    dt=tlo
    nw=0
    
    if (my_task == master_task ) then
	print *,'The total energy is ',tener0
	print *,'The total mass is ',tmass0
	print *,'The main part of this program has started'
	print *,'Number of integration steps is',iwr
    end if 
    
    
    do iws=1,iwr
	if (iws.eq.iwr) then
	    dt=tlp
	    tt=tt+tlp
	else
	    tt=tt+tlo
	end if
	nw=nw+1
	
	!if (my_task.eq. 0 .and. nt.eq.23) then
	    print *,'-------------------------------------------------------------------------------'
	    print *,'       The Energy           The Total-Mass        The iteration number'
	    nt=1
	!end if
	!------------------------------------------------
	!       The time integration
	!------------------------------------------------
	call euler(dt,iter)

	do j=1,nloc_y
	    do i=1,nloc_x
		ai=dsqrt(wh(i,j))
		u(i,j)=wu(i,j)/ai
		v(i,j)=wv(i,j)/ai
	    end do
	end do
!	if ((nyn.eq.1).or.(int(tt/43200)*43200.eq.tt)) then

	    tener=inner(wu,wv,wh,wu,wv,wh)

	    tmass = mass(wh)

!	    if (my_task.eq.0 .and. nyn.eq.0) then
		print *,'       (The integral time is      ',tt,')'
		nt=nt+1
		print *,tener,tmass,iter
		nt=nt+1
!	    endif
!	end if
    end do


    !!!!!!!!!!!!!!!!!!!!!!! Finallize  !!!!!!!!!!!!!!!!!!!!!!
    call  destroy_array
    call  MPI_BARRIER(comm, ierr)
    
end program PSWE




function inner(u1,v1,h1,u2,v2,h2)
!
    use module_para
    use module_array, only: c2
    use distribution
    use global_reductions

    implicit none

    real(r8), dimension(0:nloc_x+1,0:nloc_y+1) :: u1,v1,h1
    real(r8), dimension(0:nloc_x+1,0:nloc_y+1) :: u2,v2,h2
    real(r8), dimension(0:nloc_x+1,0:nloc_y+1) :: sum_u

    real(r8)  :: inner

    integer   :: i,j
    
    sum_u(:,:) = 0.0
    !@if (jproc == 0) then
    !@    write(*,*) u1(1:nloc_x,0)
    !@endif
    !@if (jproc == nproc_y -1) then
    !@    write(*,*) u1(1:nloc_x,nloc_y+1)
    !@endif

    do j=0,nloc_y+1
	do i=0,nloc_x+1
	    sum_u(i,j) = (u1(i,j)*u2(i,j)+v1(i,j)*v2(i,j)+h1(i,j)*h2(i,j))*c2(j)
	end do
    end do


    inner = global_sum(sum_u)

end function inner

function mass(wh1)

    use module_para
    use module_array,only: c2
    use distribution
    use global_reductions

    implicit none

    real(r8), dimension(0:nloc_x+1,0:nloc_y+1) :: wh1
    real(r8), dimension(0:nloc_x+1,0:nloc_y+1) :: whtemp

    real(r8)  :: mass

    integer   :: i,j
    
   
    do j=0,nloc_y+1
	do i=0,nloc_x+1
	    whtemp(i,j) = wh1(i,j)*c2(j)
	end do
    end do
    !write(*,*) 'c2(n) on ',my_task,  c2(nloc_y+1)
    mass = global_sum(whtemp)

end function mass
    




