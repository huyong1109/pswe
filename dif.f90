!
SUBROUTINE DIFUH(wu,wv,du,dh,h)
!
!   consider U, phi, compute advection and coriolios terms 
! 
	use module_para
	use distribution 
	use global_reductions
	use module_array, only : c1, c12, f1, f2
!
 	implicit none
!
	real(r8),dimension(0:nloc_x+1,0:nloc_y+1)  :: dh,h,hy
	real(r8),dimension(0:nloc_x+1,0:nloc_y+1)  :: wu,du,u
	real(r8),dimension(0:nloc_x+1,0:nloc_y+1)  :: wv,v,vv
	real(r8)                      :: hyn,hys
	real(r8)                      :: ff
	integer                     :: i,j
!
	do j=1,nloc_y
	   do i=1,nloc_x
	      u (i,j)=wu(i,j)/h(i,j)
	      v (i,j)=wv(i,j)/h(i,j)
	      vv(i,j)=v (i,j)*c1(j)
	   enddo
	enddo

!   advection term

	call advct(u,vv,wu,du)
!   
	do j=1,nloc_y
	   do i=1,nloc_x
	      ff=f1(j)+u(i,j)*f2(j)
	      du(i,j)=du(i,j)-ff*wv(i,j)
!
	      hy(i,j)=wv(i,j)*h(i,j)*c1(j)
	   end do
	end do
!
	!do i=2,np
	!   du(i,1)=0.0
	!   du(i,n)=0.0
!
	!   hy(i,1)=0.0
	!   hy(i,n)=0.0
	!enddo

	call polar_setval(du, 0.0, 0.0)
	call polar_setval(hy, 0.0, 0.0)
!   
	do j=1,nloc_y
	   do i=1,nloc_x
	      dh(i,j)=(hy(i,j+1)-hy(i,j-1))*c12(j)
	   end do
	end do
!
	hys=0.0
	hyn=0.0
		
	!do i=2,np
	!   hys=hys+hy(i,2)
	!   hyn=hyn-hy(i,n-1)
	!enddo
	call polar_sum(hy, hys, hyn)


	hys=hys*c12(0)/(np-1)  ! only used on proc id_y  ==  0
	hyn=hyn*c12(nloc_y+1)/(np-1) ! only used on proc id_y  == nproc_y
!
	!do i=2,np
	!   dh(i,1)=hys
	!   dh(i,n)=hyn
	!enddo
	call polar_setval(dh, hys, hyn)
!
	return
END SUBROUTINE DIFUH
!
SUBROUTINE DIFV(wu,wv,wh,dv,h)
!
    use module_para
    use distribution 
    use module_array, only : c1, c12, c14, f1, f2
!
 	implicit none
!
	real*8,dimension(0:nloc_x+1,0:nloc_y+1)  :: wh,h,hh
	real*8,dimension(0:nloc_x+1,0:nloc_y+1)  :: wu,u
	real*8,dimension(0:nloc_x+1,0:nloc_y+1)  :: wv,dv,v,vv
	real*8                      :: ff
	integer                     :: i,j
!
	do j=1,nloc_y
	   do i=1,nloc_x
	       hh(i,j)=h(i,j)*c1(j)
	       u(i,j)=wu(i,j)/h(i,j)
	       v(i,j)=wv(i,j)/h(i,j)
	       vv(i,j)=v(i,j)*c1(j)
	  enddo
	enddo
!
	call advct(u,vv,wv,dv)
!   
	do j=1,nloc_y
	   do i=1,nloc_x
	      ff=f1(j)+u(i,j)*f2(j)
	      dv(i,j)=dv(i,j)+hh(i,j)*(wh(i,j+1)-wh(i,j-1))*c12(j)+ff*wu(i,j)
	   end do
	end do
!
	!do i=2,np
	!   dv(i,1)=0.0
	!   dv(i,n)=0.0
	!enddo
	call polar_setval(dv, 0.0, 0.0)
!
	return
END SUBROUTINE DIFV
!
SUBROUTINE advct(u,v,f,df)
!
    use module_para
    use boundary
    use distribution
    use module_array, only : c1, c13, c14, f2
!
 	implicit none
!
	real*8,dimension(0:nloc_x+1,0:nloc_y+1) :: f,df
	real*8,dimension(0:nloc_x+1,0:nloc_y+1) :: u,v
	real*8,dimension(0:nloc_x+1,0:nloc_y+1) :: su,sv
	real*8                     :: dx,dy
	integer                    :: i,j
!
	do j=1,nloc_y
	   do i=1,nloc_x
	      su(i,j)=f(i,j)*u(i,j)
	      sv(i,j)=f(i,j)*v(i,j)
	   enddo
	enddo
	
	! boundary updates
	call update_ghost_cells(su)
    	call update_ghost_cells(sv)
!
	!if (jproc == 0 ) then 
	!    do i=1,nloc_y
	!       sv(i,0)=0.0
	!       df(i,0)=0.0
	!    enddo
	!endif 
	!if (jproc == nproc_y-1 ) then 
	!    do i=1,nloc_y
	!       sv(i,nloc_y+1)=0.0
	!       df(i,nloc_y+1)=0.0
	!    enddo
	!endif 
	call polar_setval(sv, 0.0, 0.0)
	call polar_setval(df, 0.0, 0.0)
	
	do j=1,nloc_y
	   do i=1,nloc_x
	      dx=u(i,j)*(f(i+1,j)-f(i-1,j))+su(i+1,j)-su(i-1,j)
	      dy=v(i,j)*(f(i,j+1)-f(i,j-1))+sv(i,j+1)-sv(i,j-1)
	      df(i,j)=dx*c13(j)+dy*c14(j)
	   enddo
	enddo
!
	return
END SUBROUTINE advct
 
subroutine polar_setval(wu, vals, valn)
    use module_para
    use distribution

    real(r8), dimension(0:nloc_x+1,0:nloc_y+1) :: wu 
    real(r8) :: vals, valn

    ! local 
    integer :: i, j 

    if(jproc == 0) then 
	do i=1,nloc_x
	  wu(i,0)=vals
 	enddo
    end if 
    if(jproc == nproc_y-1) then
	do i=1,nloc_x
	  wu(i,nloc_y+1)=valn
 	enddo
    end if 

end subroutine polar_setval
