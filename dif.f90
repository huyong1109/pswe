 subroutine DIFUH(wu,wv,du,dh,h)
!
!   consider U, phi, compute advection and coriolios terms 
!
    use module_para
    use communicate, only: master_print_message
    use module_io
    use distribution 
    use global_reductions
    use boundary
    use module_array, only : c1, c12, f1, f2

    implicit none

    real(r8),dimension(0:nloc_x+1,0:nloc_y+1)  :: dh,h,hy
    real(r8),dimension(0:nloc_x+1,0:nloc_y+1)  :: wu,du,u
    real(r8),dimension(0:nloc_x+1,0:nloc_y+1)  :: wv,v,vv
    real(r8)   :: hyn,hys
    real(r8)   :: vals,valn
    real(r8)   :: ff
    integer    :: i,j,ierr

    do j=0,nloc_y+1
        do i=0,nloc_x+1
    	u (i,j)=wu(i,j)/h(i,j)
    	v (i,j)=wv(i,j)/h(i,j)
    	vv(i,j)=v (i,j)*c1(j)
        enddo
    enddo

    
    !   advection term
    call advct(u,vv,wu,du)
    hy(:,:) = 0.0 
    do j=0,nloc_y+1
        do i=0,nloc_x+1
    	ff=f1(j)+u(i,j)*f2(j)
    	du(i,j)=du(i,j)-ff*wv(i,j)
    	
    	hy(i,j)=wv(i,j)*h(i,j)*c1(j)
        end do
    end do
    vals =0.0 
    valn =0.0 
    call polar_setval(du,vals, valn )
    call polar_setval(hy,vals, valn )
    dh(:,:) = 0.0 
    do j=1,nloc_y
        do i=1,nloc_x
    	dh(i,j)=(hy(i,j+1)-hy(i,j-1))*c12(j)
        end do
    end do
    !
    hys=0.0
    hyn=0.0
    
    call polar_sum(hy, hys, hyn)
    
    
    hys=hys*c12(0)/(np-1)  ! only used on proc id_y  ==  0
    hyn=-hyn*c12(nloc_y+1)/(np-1) ! only used on proc id_y  == nproc_y
    
    call polar_setval(dh, hys, hyn)
    call update_boundary(du)
    call update_boundary(dh)
    !
    return
 end subroutine DIFUH
!
 subroutine DIFV(wu,wv,wh,dv,h)

    use module_para
    use distribution 
    use boundary
    use module_array, only : c1, c12, c14, f1, f2

    implicit none

    real(r8),dimension(0:nloc_x+1,0:nloc_y+1), intent(in)  :: wu, wv, wh, h
    real(r8),dimension(0:nloc_x+1,0:nloc_y+1), intent(out):: dv
    real(r8),dimension(0:nloc_x+1,0:nloc_y+1)  :: hh, u, v, vv

    real(r8)                      :: ff, vals, valn
    integer                     :: i,j
    
    do j=0,nloc_y+1
        do i=0,nloc_x+1
    	hh(i,j)=h(i,j)*c1(j)
    	u(i,j)=wu(i,j)/h(i,j)
    	v(i,j)=wv(i,j)/h(i,j)
    	vv(i,j)=v(i,j)*c1(j)
        enddo
    enddo
    !
    call advct(u,vv,wv,dv)
    
    do j=1,nloc_y
        do i=1,nloc_x
    	ff=f1(j)+u(i,j)*f2(j)
    	dv(i,j)=dv(i,j)+hh(i,j)*(wh(i,j+1)-wh(i,j-1))*c12(j)+ff*wu(i,j)
        end do
    end do
    !
    vals = 0.0 
    valn = 0.0 
    call polar_setval(dv, vals, valn)
    call update_boundary(dv)
    !
    return
end subroutine DIFV
!
 subroutine advct(u,v,f,df)

    use kinds_mod, only: r8
    use communicate, only: comm
    use module_io, only: check_nan
    use boundary, only: update_boundary
    use distribution, only: nloc_x, nloc_y
    use module_array, only : c1, c13, c14, f2
    !
    implicit none
    !
    real(r8),dimension(0:nloc_x+1,0:nloc_y+1), intent(inout):: df
    real(r8),dimension(0:nloc_x+1,0:nloc_y+1), intent(in) :: u,v,f 
    real(r8),dimension(0:nloc_x+1,0:nloc_y+1) :: su,sv
    real(r8)                     :: dx,dy
    real(r8)                     :: vals, valn
    integer                    :: i,j,ierr
    !
    do j=0,nloc_y+1
        do i=0,nloc_x+1
    	su(i,j)=f(i,j)*u(i,j)
    	sv(i,j)=f(i,j)*v(i,j)
        enddo
    enddo
    df(:,:) = 0.0 
    
    vals = 0.0
    valn = 0.0
    call polar_setval(sv, vals, valn)
    call polar_setval(df, vals, valn)
    
    do j=1,nloc_y
        do i=1,nloc_x
    	dx=u(i,j)*(f(i+1,j)-f(i-1,j))+su(i+1,j)-su(i-1,j)
    	dy=v(i,j)*(f(i,j+1)-f(i,j-1))+sv(i,j+1)-sv(i,j-1)

    	df(i,j)=dx*c13(j)+dy*c14(j)
        enddo
    enddo

    call update_boundary(df)
    call MPI_BARRIER(comm, ierr)
    return
 end subroutine advct

 subroutine polar_setval(wu, vals, valn)
    use kinds_mod, only: r8
    use module_para, only: nproc_y
    use module_io, only: check_nan
    use distribution, only: jproc, nloc_x, nloc_y
    
    real(r8), dimension(0:nloc_x+1,0:nloc_y+1), intent(inout) :: wu 
    real(r8), intent(in):: vals, valn
    
    ! local 
    integer :: i, j 
    
    if(jproc == 0) then 
        do i=1,nloc_x
	    wu(i,0)=vals
        enddo
    endif 

    if(jproc == nproc_y-1) then
        do i=1,nloc_x
	    wu(i,nloc_y+1)=valn
        enddo
    endif 
    call check_nan(wu, 'wu in polar_setval')
    return
 end subroutine polar_setval
