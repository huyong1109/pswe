! ---------------------------------------------------------------------------------------!
! This subroutine is to provide initial condition using four-wave Rossby-Haurwitz waves	 !
! ---------------------------------------------------------------------------------------!
!
SUBROUTINE haurwitz
!
    use module_para
    use module_array
    use distribution
    use boundary
!
    implicit none
!
! 	real*8,  parameter ::    omg  = 7.848d-6		 ! angular velocity of RH wave
 	real*8,  parameter ::    omg  = 3.924d-6		 ! angular velocity of RH wave
    real*8,  parameter ::    fi0  = 78400d0			 ! minimun potential height
	real*8,  parameter ::    r    = 4d0              ! wave number of RH wave
!													 !
    real*8             ::   af,ai,aj,ak,al			 ! working variable
    real*8             ::   bf,bi,bj,fi,r1,r2		 ! working variable
    real*8             ::   cf,detar,u0,v0,u1   	 ! working variable
!													 !
    integer            ::   i,j                      ! working variable
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
	      bi=i*detar-detar
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
!
	do j = 1, nloc_y
	   wh(1,j)=wh(np,j)
	   wh(np+1,j)=wh(2,j)
!
	   u(1,j)= u(np,j)
	   u(np+1,j)= u(2,j)
!
	   v(1,j)= v(np,j)
	   v(np+1,j)= v(2,j)
	end do
!
    if(jproc == 0) then 
	do i=1,nloc_x
	  fi = fi0
	  wh(i,0)=fi
	  u(i,0)=0
	  v(i,0)=0
 	enddo
    end if 
    if(jproc == nproc_y-1) then
	do i=1,nloc_x
	  fi = fi0
	  wh(i,nloc_y+1)=fi
	  u(i,nloc_y+1)=0
	  v(i,nloc_y+1)=0
 	enddo
    end if 

    call update_ghost_cells(wh)
    call update_ghost_cells(u)
    call update_ghost_cells(v)
!
	return
END SUBROUTINE haurwitz
