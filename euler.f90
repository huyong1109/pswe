!
subroutine euler(dt,iter)

    use module_para
    use communicate, only:master_print_message
    use module_array, only: wu, wv, wh, c11
    use module_io
    use boundary, only: update_boundary
    use global_reductions, only: global_sum

    implicit none

    integer i,i1,i2,j,k,iter, ierr, break_cir

    real(r8),dimension(0:nloc_x+1,0:nloc_y+1)  :: tu,tv,th,h 
    real(r8),dimension(0:nloc_x+1,0:nloc_y+1)  :: du,dv,dh 

    real(r8),dimension(0:nloc_x+1,0:nloc_y+1)      :: a1,a2,rh,ru 

    real(r8)                      :: ai,aj,dt,dt2,en,en0,den

    real(r8),external             :: inner

    !	half timestep
    dt2=dt*0.5d0
    call update_boundary(wu)
    call update_boundary(wv)
    call update_boundary(wh)
    do j=0,nloc_y+1
	do i=0,nloc_x+1
	    tu(i,j)=wu(i,j)
	    tv(i,j)=wv(i,j)
	    th(i,j)=wh(i,j)
	end do
    end do

    !	Initial energy
    en0=inner(wu,wv,wh,wu,wv,wh)

    !	Iterate until convergence ( |E-E0| < 1.0E-15 E0 )
    break_cir = 0
    do k=1,1000

	do j=0,nloc_y+1
	    do i=0,nloc_x+1
		h(i,j)=dsqrt(th(i,j))
	    end do
	end do
	
	!	consider U, phi, compute advection and coriolios terms 

	call DIFUH(tu,tv,du,dh,h)
	!	get the tri-diagonal matrix
	do j=0,nloc_y+1
	    do i=0,nloc_x+1
		tu(i,j)=wu(i,j)-dt2*du(i,j)
		th(i,j)=wh(i,j)-dt2*dh(i,j)
	    end do
	end do
	
	a1(:,:) = 0.0
	a2(:,:) = 0.0
	do j=1,nloc_y
	    ai=dt2*c11(j)

	    do i=0,nloc_x+1
		aj=ai*h(i,j)
		a1(i,j)=aj
		ru(i,j)=tu(i,j)*aj
		a2(i,j)=aj*aj
	    end do
	end do 
	call update_boundary(ru) 
		
	do j=0,nloc_y+1
	    do i=0,nloc_x+1
		rh(i, j)=th(i,j)-ru(i+1, j)+ru(i-1,j )
	    enddo
	enddo
	call  Tri_diag(a2, rh, th)

	call update_boundary(th) 

	do j = 1, nloc_y
	    do i = 1, nloc_x
		tu(i,j)=tu(i,j)-a1(i,j)*(th(i+1,j)-th(i-1,j))
	    enddo
	enddo 

	call update_boundary(tu) 
	
	call DIFV(tu,tv,th,dv,h)
	
	do j=1,nloc_y
	    do i=1, nloc_x
		tv(i,j)=wv(i,j)-dt2*dv(i,j)
	    enddo
	enddo
	
	call update_boundary(tv) 
	
	en=inner(tu,tv,th,tu,tv,th)
	
	den=dabs(en-en0)*2.0/(en+en0)
	!if (my_task .eq. master_task ) write(*,*) k ,  en
	en0=en
	if (den.lt.1.0d-15) break_cir = 1
	if (break_cir .eq. 1) goto 10
    enddo

    10   continue
    iter=k


    do j=0,nloc_y+1
	do i=0,nloc_x+1
	    wu(i,j)=tu(i,j)*2.0d0-wu(i,j)
	    wv(i,j)=tv(i,j)*2.0d0-wv(i,j)
	    wh(i,j)=th(i,j)*2.0d0-wh(i,j)
	end do
    enddo

 end subroutine euler

 subroutine Tri_diag(ra, rh, th )
    use kinds_mod, only: r8
    use distribution, only: nloc_x, nloc_y, jproc, iproc
    use module_para, only: nproc_x, kn
    use module_io, only: check_nan
    use boundary, only: update_boundary
    use global_reductions, only:  lat_gather, lat_scatter


    implicit none

    integer i,i1,i2,j,k,iter,ierr

    real(r8),dimension(0:nloc_x+1,0:nloc_y+1), intent(in)    :: ra,rh
    real(r8),dimension(0:nloc_x+1,0:nloc_y+1), intent(inout)   :: th
    real(r8),dimension(:,:), allocatable        :: fm,fp,f0,gm,gp,g0,rf,rg 
    real(r8),dimension(:,:), allocatable	      ::  gra, grh, gth

    real(r8)                      :: ai,aj


    if (iproc == 0 ) then
        allocate(gth(1:nloc_x, 1:nloc_y*nproc_x))
        allocate(gra(1:nloc_x, 1:nloc_y*nproc_x))
        allocate(grh(1:nloc_x, 1:nloc_y*nproc_x))
       
	allocate(fm(1:kn,1:nloc_y))
        allocate(fp(1:kn,1:nloc_y))

        allocate(f0(1:kn,1:nloc_y))
        allocate(gm(1:kn,1:nloc_y))
        allocate(gp(1:kn,1:nloc_y))
        allocate(g0(1:kn,1:nloc_y))
        allocate(rf(1:kn,1:nloc_y))
        allocate(rg(1:kn,1:nloc_y))


        gth(:,:) = 0.0
        gra(:,:) = 0.0
        grh(:,:) = 0.0
    endif 

    call lat_gather(gra, ra)
    call lat_gather(grh, rh)

    !	grid divided into two groups. 
    if (iproc == 0 ) then

        do j = 1, nloc_y
            do i=1,kn
                i1=i*2-1
                i2=i1+1

                fp(i,j)=-gra(mod(i2-1, nloc_x)+1, j+i1/nloc_x*nloc_y)
                rf(i,j)=grh(mod(i1-1, nloc_x)+1,j+i1/nloc_x*nloc_y)

                gm(i,j)=-gra(mod(i1-1, nloc_x)+1,j+i1/nloc_x*nloc_y)
                rg(i,j)=grh(mod(i2-1, nloc_x)+1, j+i1/nloc_x*nloc_y)
            enddo

	    do i=2,kn
		fm(i,j)=fp(i-1,j)
	    enddo
	    fm(1,j)=fp(kn,j)

	    do i=1,kn-1
		gp(i,j)=gm(i+1,j)
	    enddo
	    gp(kn,j)=gm(1,j)

	    do i=1,kn
		f0(i,j)=1.0-fm(i,j)-fp(i,j)
		g0(i,j)=1.0-gm(i,j)-gp(i,j)
	    enddo
	enddo


	    call LU(fm(:,:),f0(:,:),fp(:,:),rf(:,:),nloc_y, kn)
	    call LU(gm(:,:),g0(:,:),gp(:,:),rg(:,:),nloc_y, kn)


	do j=1,nloc_y
	    do i=1,kn
		i1=i*2-1
		i2=i1+1

		gth(mod(i1-1, nloc_x)+1,j+i1/nloc_x*nloc_y)=rf(i,j)
		gth(mod(i2-1, nloc_x)+1, j+i1/nloc_x*nloc_y)=rg(i,j)
	    enddo
	enddo


    endif

    call lat_scatter(gth, th)	
    if (iproc == 0 ) then
        deallocate(fm)
        deallocate(fp)
        deallocate(f0)
        deallocate(gm)
        deallocate(gp)
        deallocate(g0)
        deallocate(rf)
        deallocate(rg)
        deallocate(gth)
        deallocate(gra)
        deallocate(grh)
    endif 
 
 end subroutine Tri_diag
 subroutine LU(a,b,c,r,m,n)
    !Gauss elimination method/Chasing method  with periodic boundary
    ! on Wang's book P159
    use kinds_mod
    implicit none
    !
    integer                 :: i,j, n, m
    real(r8) ,dimension(1:m)     :: sn,rn
    real(r8)		         :: ai
    !
    real(r8),dimension(1:n, 1:m) :: a,b,c,r
    real(r8),dimension(1:n, 1:m) :: s,t
    !
    s(1, 1:m)=a(1, 1:m)  ! s = \tilde(a)
    t(1, 1:m)=c(n, 1:m)  ! t = \tilde(c)
    !
    sn(:)=0.0	   ! \tilde(b)_N
    rn(:)=0.0	   ! \tilde(r)_N
    
    do j = 1, m
	do i=2,n-1
	    ai=a(i, j)/b(i-1, j)
	    b(i, j)=b(i, j)-ai*c(i-1,j)
	    r(i, j)=r(i, j)-ai*r(i-1,j)
	    s(i, j)=-ai*s(i-1, j)
	    !
	    ai=t(i-1, j)/b(i-1, j)
	    t(i, j)=-ai*c(i-1, j)
	    sn(j)  =sn(j)-ai*s(i-1, j)
	    rn(j)  =rn(j)-ai*r(i-1, j)
	enddo
    enddo
    
    
    do j = 1, m
	a(n, j)=a(n, j)+t(n-1, j)
    	b(n, j)=b(n, j)+sn(j)
    	c(n-1, j)=c(n-1, j)+s(n-1, j)
    	r(n, j)=r(n, j)+rn(j)

	ai=a(n, j)/b(n-1, j)
    	b(n, j)=b(n, j)-ai*c(n-1, j)
    	r(n, j)=r(n, j)-ai*r(n-1, j)
    	!
    	r(n, j)=r(n, j)/b(n, j)
    	r(n-1, j)=(r(n-1, j)-c(n-1, j)*r(n, j))/b(n-1, j)

	do i=n-2,1,-1
	    ai=r(i, j)-s(i, j)*r(n, j)-c(i, j)*r(i+1, j)
	    r(i, j)=ai/b(i, j)
	enddo
    enddo 

    return
    
 end subroutine LU
    
 subroutine LU0(a,b,c,r,n)
    !Gauss elimination method/Chasing method  with periodic boundary
    ! on Wang's book P159
    use kinds_mod
    implicit none
    !
    integer                 :: i,n
    real(r8)                :: ai,sn,rn
    !
    real(r8),dimension(1:n) :: a,b,c,r
    real(r8),dimension(1:n) :: s,t
    !
    s(1)=a(1)  ! s = \tilde(a)
    t(1)=c(n)  ! t = \tilde(c)
    !
    sn=0.0	   ! \tilde(b)_N
    rn=0.0	   ! \tilde(r)_N
    !
    do i=2,n-1
	ai=a(i)/b(i-1)
	b(i)=b(i)-ai*c(i-1)
	r(i)=r(i)-ai*r(i-1)
	s(i)=-ai*s(i-1)
	!
	ai=t(i-1)/b(i-1)
	t(i)=-ai*c(i-1)
	sn  =sn-ai*s(i-1)
	rn  =rn-ai*r(i-1)
    enddo
    !
    a(n)=a(n)+t(n-1)
    b(n)=b(n)+sn
    c(n-1)=c(n-1)+s(n-1)
    r(n)=r(n)+rn
    !
    ai=a(n)/b(n-1)
    b(n)=b(n)-ai*c(n-1)
    r(n)=r(n)-ai*r(n-1)
    !
    r(n)=r(n)/b(n)
    r(n-1)=(r(n-1)-c(n-1)*r(n))/b(n-1)
    !
    do i=n-2,1,-1
	ai=r(i)-s(i)*r(n)-c(i)*r(i+1)
	r(i)=ai/b(i)
    enddo
    !
    return
    !
 end subroutine LU0


