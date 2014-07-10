! ---------------------------------------------------------------------------------------!
! This subroutine is to provide initial condition using four-wave Rossby-Haurwitz waves	 !
! ---------------------------------------------------------------------------------------!
!
 subroutine haurwitz

    use module_para
    use communicate, only: master_print_message
    use module_array
    use distribution
    use boundary
    use module_io

    implicit none

    !real*8,  parameter ::    omg  = 7.848d-6		 ! angular velocity of RH wave
    real*8,  parameter ::    omg  = 3.924d-6		 ! angular velocity of RH wave
    real*8,  parameter ::    fi0  = 78400d0			 ! minimun potential height
    real*8,  parameter ::    r    = 4d0              ! wave number of RH wave
    !													 !
    real*8             ::   af,ai,aj,ak,al			 ! working variable
    real*8             ::   bf,bi,bj,fi,r1,r2		 ! working variable
    real*8             ::   cf,detar,u0,v0,u1   	 ! working variable
    real*8             ::   vals, valn
    !													 !
    integer            ::   i,j, ierr                     ! working variable
    !
    detar=2*pi/p*r
    r1=r+1
    r2=r*r
    !
    do j=1,nloc_y
        do i=1,nloc_x
    	aj=c1(j)
    	ai=s1(j)
    	!-------------------------- U(x,y,0) ------------------------------
    	ak=aj**r
    	al=aj*aj
    	bi=(i+iglobal+1)*detar-detar
    	bj=dcos(bi)
    	u1=aj+ak/aj*ai*ai*bj*r-ak*aj*bj
    	u0=u1*a*omg
    	!--------------------------- V(x,y,0) -----------------------------
    	bj=dsin(bi)
    	v0=-a*r*omg*ak/aj*ai*bj
    	!--------------------------- H(x,y,0) -----------------------------
    	bj=dcos(bi*2)
    	bi=dcos(bi)
    	af=r1*al+2*r2-r-2-2*r2/al
    	af=af*ak*ak
    	af=af*omg*omg/4+omg*(omg0*2+omg)*al/2
    	bf=r2+2*r+2-r1*r1*al
    	bf=bf*ak*2*omg*(omg+omg0)
    	bf=bf/r1/(r+2)
    	cf=r1*al-r-2
    	cf=cf*omg*omg*ak*ak/4
    	fi=af+bf*bi+cf*bj
    	fi=fi0+fi*a*a
    	!------------------------------------------------------------------
    	wh(i,j)=fi
    	u(i,j)=u0
    	v(i,j)=v0
        enddo
    enddo
    call polar_setval(wh, fi0, fi0)
    vals = 0.0
    valn = 0.0
    call polar_setval(u,vals,valn )
    call polar_setval(v,vals,valn )
    
 
    
    
    call update_boundary(wh)
    call update_boundary(u)
    call update_boundary(v)
    !write(*,*) wh
    !call mpi_print(v)
    !
    call master_print_message('Finish Haurwitz')
    return
 end subroutine haurwitz
